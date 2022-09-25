function out_data = alg_radp_matched(alg_settings)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% ROBUST ADAPTIVE DYNAMIC PROGRAMMING (RADP) ALGORITHM -- MATCHED
% UNCERTAINTY
%
% Brent Wallace  
%
% 2022-01-17
%
% This program implements the RADP algorithm (with matched uncertainty)
% presented in,
%
%   Y. Jiang and Z.-P. Jiang. "Robust adaptive dynamic programming and
%   feedback stabilization of nonlinear systems." IEEE Transactions on
%   Neural Networks and Learning Systems, 25:882-893, 2014.
%
% *** NOTE:
%
% This program was developed to run single-input systems (m = 1) ONLY. It
% can easily be adapted to the multi-input case (cf. Wallace, Si TNNLS 2022
% Part I, Sec. VI).
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% out_data = alg_radp_matched(alg_settings)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% alg_settings  struct with the following fields:
%   
%   preset                  (String) example preset (see main.m for
%                           options).
%   sys                     (Struct) contains system tag/info. See notes in
%                           'config.m' for specific fields.
%   alg                     (String) redundant for this function. Contains
%                           the tag of this algorithm.
%   Q                       (n x n matrix, or string) If matrix, is the
%                           positive definite state penalty matrix. If
%                           string, is the tag of the desired non-quadratic
%                           positive definite state penalty function.
%   R                       (m x m matrix, or string) If matrix, is the
%                           positive definite control penalty matrix. If
%                           string, is the tag of the desired non-quadratic
%                           positive definite control penalty function.
%   basis                   (Struct) contains activation function bases
%                           parameters. Has the following fields:
%       .Phi               (Struct) contains critic NN activation function
%                           parameters (cf. eqn. (13)). Has fields:
%           .tag            (String) tag of desired activation functions.
%           .N              (Integer) number of critic NN activation
%                           functions N1.
%       .Psi               (Struct) contains actor NN activation function
%                           parameters (cf. eqn. (14)). Has fields:
%           .tag            (String) tag of desired activation functions.
%           .N              (Integer) number of actor NN activation
%                           functions N2.
%   noise                   (Struct) contains info for probing noise. Has
%                           the following fields:
%       .tag                (String) tag of specific probing noise signal
%                           to be injected (see eval_noise.m for options).
%   do_uncertainty          (Bool) 1 = include dynamic uncertainty (i.e.,
%                           the w-subsystem, cf. eqn. (8)). 0 = don't
%                           include dynamic uncertainty.
%   rho                     (String, optional) tag corresponding to the
%                           robust redesign function \rho(s) (cf. under
%                           eqn. (31)). NOTE: If system does not have
%                           uncertainty, then declare this field empty.
%   l                       (Integer) number of samples to collect for each
%                           least-squares minimization (cf. eqn. (15)).
%   tf                      (Double) Length of learning window [0, t_f]
%                           (sec).
%   tsim                    (Double) Length of simulation to run after
%                           learning window (sec). I.e., post-learning
%                           simulation happens [t_f, t_f + tsim].
%   u_0                     (String) Tag corresponding to the initial
%                           stabilizing policy u_0(x) (cf. Assumption 3.2).
%                           Evaluation of this policy is handled by
%                           eval_u.m, so the policy needs to be encoded
%                           there.
%   w0                      (p-dimensional vector, optional) ICs for
%                           unknown dynamics w(t_0). NOTE: If system does
%                           not have uncertainty, then declare this field
%                           empty.
%   x0                      (n-dimensional vector) ICs for known x-dynamics
%                           x(t_0)..
%   w_0                     (N2-dimensional vector) ICs for actor NN
%                           weights.
%   eps1                    (Double, optional) Tolerance \epsilon_1 to
%                           determine phase-1 learning termination (cf.
%                           eqn. (59)). NOTE: This field need only be
%                           declared if the termination condition selected
%                           is convergence-based (i.e., if term_mode ==
%                           'conv', see below).
%   istar                   (Integer, optional) Number of iterations to run
%                           algorithm for before manually terminating.
%                           NOTE: This field need only be declared if the
%                           termination condition selected is
%                           manual termination (i.e., if term_mode ==
%                           'manual', see below).
%   term_mode               (String) Termination mode of the algorithm. Has
%                           the following options:
%       'conv'              Convergence-based termination (as determined by
%                           threshold 'eps1').
%       'manual'            Manual termination after iteration count
%                           'istar'.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% out_data                  (Struct) algorithm output data. Has the
%                           following fields:
%       .tvec               ('simlength'-dimensional vector) vector of
%                           'simlength' time indices corresponding to the
%                           simulation time instants over the course of the
%                           algorithm execution.
%       .xmat               ('simlength' x (p+n+1) matrix) Matrix whose row
%                           indexes the time instants specified in .tvec,
%                           and whose n-columns are the state vector at the
%                           respective time instant.
%           xmat(k,:) =   [     w(tvec(k))   R^p
%                               x(tvec(k))   R^n   ]^T
%       .umat               ('simlength' x 1 matrix) Matrix (vector,
%                           really) whose row indexes the time instants
%                           specified in .tvec, and whose 1 column is the
%                           control signal u(t) at the respective time
%                           instant.
%       .istar              (Integer) final termination index of the PI
%                           algorithm.
%       .c_mat            ('istar'+1 x N1 matrix) Matrix whose row
%                           indexes the PI index, and whose N1-columns are
%                           the critic NN weights c_i at the respective PI
%                           index i.
%       .w_mat            ('istar'+1 x N2 matrix) Matrix whose row
%                           indexes the PI index, and whose N2-columns are
%                           the actor NN weights w_i at the respective PI
%                           index i.
%       cond_A_vec          ('istar'-dim. vector) The i-th index of
%                           this vector contains the condition number of
%                           the matrix involved in performing the
%                           least-squares minimization associated with the
%                           i-th iteration weight update.
%       
%
% *************************************************************************
% *************************************************************************
% *************************************************************************



% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INITIALIZE
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
% 
% GLOBAL VARIABLES
% 
% *************************************************************************


global preset;
global sys;

global Q;
global R;

global basis;
global noise;

global u_0;
% global l;
% global T;

% Keeps track of if the algorithm is in the learning phase or the final
% phase
global is_learning;

global rho;

% Keeps track of if dynamic uncertainty is to be included in this problem.
global do_uncertainty;

% Critic NN parameters
global c_i;

% Actor NN parameters
global w_i;




% *************************************************************************
% 
% UNPACK ALGORITHM SETTINGS/PARAMETERS
% 
% *************************************************************************


preset = alg_settings.preset;               
sys = alg_settings.sys;         % System

do_uncertainty = alg_settings.do_uncertainty;

% System dimensions
n = alg_settings.sys.n;         % Order of known x-dynamics
if do_uncertainty
    p = alg_settings.sys.p;         % Order of unknown w-dynamics
end

Q = alg_settings.Q;
R = alg_settings.R;

noise = alg_settings.noise;

u_0 = alg_settings.u_0;

rho = alg_settings.rho;

l = alg_settings.l;

% Simulation times
tf = alg_settings.tf;
tsim = alg_settings.tsim;

% Algorithm termination mode
term_mode = alg_settings.term_mode;

% Extract termination parameters
switch term_mode

    case 'conv'

        eps1 = alg_settings.eps1;

    case 'manual'

        istar = alg_settings.istar;

    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    otherwise

        error('**** ERROR: TERMINATION MODE TAG NOT RECOGNIZED ***');              
end

% ***********************
%       
% BASIS SETTINGS
%
% See description for details.
%

% Basis struct
basis = alg_settings.basis;

% Critic basis for approximating cost function V_i(x) (cf. eqn. (13))
N1 = alg_settings.basis.Phi.N;    
% Actor basis for approximating policy u_{i + 1}(x) (cf. eqn. (14))
N2 = alg_settings.basis.Psi.N;



% *************************************************************************
% 
% ALGORITHM INITIALIZATION
% 
% *************************************************************************

% ***********************
%       
% INITIAL CONDITIONS
%
% NOTE: State partitioning for simulation
%
% In order to perform the RADP algorithm with unmatched uncertainty, the
% following three types of dynamic variables will be necessary for
% simulation:
%
%   System states (including w, x, cf. eqns. (8)-(9))
%   Integrals associated with the least squares minimization (15)
%

% Initial conditions
w0 = alg_settings.w0;
xs0 = alg_settings.x0;


% Initial condition for simulation. See odefunct() for a description of the
% state partition.
if do_uncertainty
    x0_sim = [  w0
                xs0
                zeros(N2^2, 1)
                zeros(N2, 1);
                0               ];
else
    x0_sim = [  xs0
                zeros(N2^2, 1)
                zeros(N2, 1);
                0               ];    
end


% ***********************
%       
% MISCELLANEOUS VARIABLES
%

% Initialize critic NN weight vector
c_i = zeros(N1, 1);     % Holds critic NN weights of current iteration

% Initialize actor NN weight vector
w_i = alg_settings.w_0; % Holds actor NN weights of current iteration   

% Initialize vector to store previous iteration critic NN weights
c_im1 = zeros(N1, 1);   % Holds critic NN weights of previous iteration
c_im1(1) = Inf;         % So PI algorithm starts properly

% % Initialize vector to store previous iteration actor NN weights
% w_im1 = zeros(N2, 1);   % Holds critic NN weights of previous iteration
% w_im1(1) = Inf;         % So PI algorithm starts properly

% Sample collection window T
T = tf / l;

% Set learning flag
is_learning = 1;

% Initialize PI index
i = 1;


% Initialize vector to hold condition number at each iteration
cond_A_vec = [];

% *************************************************************************
% 
% DATA STORAGE
% 
% *************************************************************************

% ***********************
%       
% LEAST SQUARES MINIMIZATION VARIABLES
%
% The dynamic variables
%
%   IPsiPsi(k)   R^{N_2^2}
%   IPsiupd(k)    R^{N_2}
%   IQ(k)          R^1      
%
%       0 <= k <= l-1
%
% (see odefunct() for a description of these variables) need to be stored
% in corresponding l-row matrices (or l-dimensional vectors) for subsequent
% least squares minimization at the end of the current PI iteration. These
% storage variables are initialized as zeros here.
%

IPsiPsi_mat = zeros(l,N2^2);
IPsiupd_mat = zeros(l,N2);
IQ_vec = zeros(l,1);


% In addition, the non-dynamic variables
%
%   d_Phi(k) = \Phi_{1}(x(t_{k+1})) - \Phi_{1}(x(t_k))       R^{N_1}
%
%       0 <= k <= l-1
%
% (see the description of odefunct() for definition of the function
% \Phi_{1}) need to be stored.
%

d_Phi_mat = zeros(l, N1);


% ***********************
%       
% STATE VARIABLES
%

% Time vector corresponding to the total simulation (including the PI phase
% and the final approximate optimal policy phase)
tvec = [];

% Matrix consisting of the system state vector in each row, corresponding
% to the times in tvec
xmat = [];


% ***********************
%       
% WEIGHTS
%

% Stores actor weights c_i at each iteration of the PI algorithm (iteration
% indexed by row)
c_mat = c_i';

% Stores actor weights w_i at each iteration of the PI algorithm (iteration
% indexed by row)
w_mat = w_i';

% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% BEGIN MAIN
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

%%
% *************************************************************************
% *************************************************************************
%
% PHASE 1, PHASE 2 LEARNING
%
% See steps 1-3 of Algorithm 1, Jiang and Jiang
% 
% *************************************************************************
% *************************************************************************


% *************************************************************************
%
% COLLECT SIMULATION DATA
%
% This loop collects the l data points needed to perform the least-squares
% minimization (15).
%  
% *************************************************************************

for k = 0:l-1

    % ***********************
    %       
    % RUN SIMULATION
    %
 
    tspan = k * T + [0, T];

    [t, x] = ode45(@odefunct, tspan, x0_sim);

%     % DEBUGGING: Plot state trajectory
%     figure(100)
%     hold on
%     if do_uncertainty
%         plot(t,x(:,1:p+n));
%     else
%         plot(t,x(:,1:n));
%     end
%     title('State Trajectory')
%     grid on
%     xlabel('Time (sec)')
%     ylabel('x(t)')

    % ***********************
    %       
    % EXTRACT STATE DATA
    %         

    % Extract dynamic state variable data at beginning of simulation
    x0 = (x(1,:))';

    % Extract dynamic state variable data at end of simulation
    x1 = (x(size(x, 1),:))';

    % Unpack state data at beginning of simulation. See odefunct() for
    % a description of the state partition.
    
    % State dynamic variables
    if do_uncertainty
        w0 = x0(1:p);               % Unknown w-dynamics
        xs0 = x0(p+1:p+n);          % Known x-dynamics
    else
        xs0 = x0(1:n);              % Known x-dynamics
    end  

    % Unpack state data at end of simulation. See odefunct() for a
    % description of the state partition.
    
    % State dynamic variables
    if do_uncertainty
        w1 = x1(1:p);               % Unknown w-dynamics
        xs1 = x1(p+1:p+n);          % Known x-dynamics
    else
        xs1 = x1(1:n);              % Known x-dynamics
    end

    % Dynamic variables associated with least squares minimizations
    if do_uncertainty
        x1_lsq = x1(p+n+1:end);
    else
        x1_lsq = x1(n+1:end);
    end 

    IPsiPsi = x1_lsq(1:N2^2);
    IPsiupd = x1_lsq(N2^2+1:N2^2+N2);
    IQ = x1_lsq(N2^2+N2+1);

    % ***********************
    %       
    % DATA STORAGE -- SYSTEM STATE DATA, DYNAMIC LEAST-SQUARES DATA
    % 

    % Store time data
    tvec = [    tvec
                t       ];

    % Store system state data 
    if do_uncertainty
        xmat = [    xmat
                    x(:,1:p+n)  ];
    else
        xmat = [    xmat
                    x(:,1:n)    ]; 
    end


    % Store dynamic variables associated with least squares data 
    IPsiPsi_mat(k+1,:) = IPsiPsi';
    IPsiupd_mat(k+1,:) = IPsiupd';
    IQ_vec(k+1) = IQ;


    % ***********************
    %       
    % DATA STORAGE -- REMAINING LEAST SQUARES DATA
    %         

    % Evaluate critic basis functions \Phi_{1}(x) = (\phi_1(x), ...,
    % \phi_{N1}(x)) at x = x(t_k) and x = x(t_{k+1})
    Phix_tk = eval_phi(xs0, basis.Phi);
    Phix_tkp1 = eval_phi(xs1, basis.Phi); 

    % Calculate and store \Phi_{1}(x(t_{k+1})) - \Phi_{1}(x(t_k))
    d_Phi_mat(k+1,:) = (Phix_tkp1 - Phix_tk)';

    % ***********************
    %       
    % PREPARE FOR NEXT SIMULATION
    % 

    % IC for next simulation. System state variables for next
    % simulation are carried over from their final values in this
    % simulation. Integration variables associated with least squares
    % minimizations are reset to zero.  
    %
    % See odefunct() for a description of the state partition.
    %
    if do_uncertainty
        x0_sim = [  w1
                    xs1
                    zeros(N2^2, 1)
                    zeros(N2, 1);
                    0               ];
    else
        x0_sim = [  xs1
                    zeros(N2^2, 1)
                    zeros(N2, 1);
                    0               ];    
    end    
         

end


% *************************************************************************
%
% PHASE 1 LEARNING: PI LOOP
%
% *************************************************************************

% This bool controls when the algorithm will exit the learning loop
do_loop = 1;

while do_loop
    
    % DEBUGGING: Show iteration count, norm difference
    disp('*****')
    disp(['i = ' num2str(i) ',     ||c_i - c_{i-1}|| =     ' ...
        num2str(norm(c_i - c_im1))])
    
    % *********************************************************************
    %       
    % LEAST SQUARES MINIMIZATION (15)
    %
    % Defined as in the description of odefunct(), the least squares
    % minimization (15) has the form:
    %
    %   [A1 A2] [ c_i    = - v
    %             w_i ]
    %
    % Expressed in terms of the variables we have declared above, this is:
    %
    % A1 = d_Phi_mat;
    % A2 = 2r * ( IPsiupd_mat - IPsiPsi_mat * kron(w_i, eye(N2)) )
    % v = IQ_vec + IPsiPsi_mat * kron(w_i, w_i)
    %
    
    % Least squares variables
    A1 = d_Phi_mat;
    A2 = 2 * R * ( IPsiupd_mat - IPsiPsi_mat * kron(w_i, eye(N2)) );
    v = IQ_vec + IPsiPsi_mat * R * kron(w_i, w_i);
    
    % Perform least squares
    A = [ A1, A2 ];
    lsq = A\(-v);
    
    % Store condition number of matrix involved in least-squares
    % minimization
    cond_A_vec = [cond_A_vec ; cond(A)];
    
    % Store current w_i as w_{i-1}, extract new least squares solutions
    c_im1 = c_i;
%     w_im1 = w_i;
    c_i = lsq(1:N1);
    w_i = lsq(N1+1:end);
    
    % Store c_i, w_i for later plotting
    c_mat =   [   c_mat
                    c_i'        ];
    w_mat =   [   w_mat
                    w_i'        ];
                
    % DEBUGGING: Check problem conditioning
    disp(['Condition Number of "A" for Least Squares:           '...
                num2str(cond(A), 4)])
    disp(['Condition Number of "A^T A" for Least Squares:       '...
                num2str(cond(A'*A), 4)])
    
    % DEBUGGING: Check PE condition
    thetaTtheta_i = zeros(N1+N2);
    for ik = 1:l
        theta_ik = [A1(ik,:), A2(ik,:)];
        thetaTtheta_i = thetaTtheta_i + theta_ik' * theta_ik;
    end
    eigs = eig(thetaTtheta_i);
    disp(['Min. E-Val of \theta^T * \theta:     ' num2str(min(eigs))])
                
    % *********************************************************************
    %       
    % PREPARE FOR NEXT PI ITERATION
    %    
    
    % Check if algorithm should continue to next iteration
    switch term_mode
        
        case 'conv'
        
            do_loop = norm(c_i - c_im1) > eps1;
            
        case 'manual'
            
            do_loop = i < istar;
    
        % THROW ERROR IF TAG DOES NOT COME UP A MATCH
        otherwise
        
            error('**** ERROR: TERMINATION MODE TAG NOT RECOGNIZED ***');              
    end
    
    % Increment PI index if next iteration is to be performed
    if do_loop
        
        i = i + 1;
        
    end
    
    
end

% DEBUGGING: ICs
x_0 = alg_settings.x0

% DEBUGGING: Final critic weights c_i
c_i

% DEBUGGING: Final actor weights w_i
w_i


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% APPLY APPROXIMATE ROBUST OPTIMAL CONTROL POLICY (31)
%
% See steps 4-5 of Algorithm 1, Jiang and Jiang
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% Clear learning flag
is_learning = 0;

% ***********************
%       
% RUN SIMULATION
%

% Initial condition. Note now that the least-squares dynamic variable
% states are no longer being simulated, so the state is simply partitioned
% as, 
%
% WITH UNCERTAINTY:
%
% x_sim =   [   w           R^p
%               x           R^n     ]
%
% WITHOUT UNCERTAINTY:
%
% x_sim =   [   x           R^n     ]
%
if do_uncertainty
    x0_sim = [  w1
                xs1     ];
else
    x0_sim = xs1;   
end    
    
% Time span for simulation
tspan = [tf, tf + tsim];

% % DEBUGGING: Manual ICs
% x0_sim = [0.01 ; 0 ; 0.01 ; 0];

% Run simulation
[t, x] = ode45(@odefunct_final, tspan, x0_sim);


% ***********************
%       
% STORE DATA
%

% Store how long tvec is at this point for later
len_tvec_learn = size(tvec,1);

% Store time data
tvec = [    tvec
            t       ];


% Store system state data
xmat = [    xmat
            x       ];
            
%%            
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% PREPARE OUTPUT DATA
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% Trajectory data
out_data.tvec = tvec;
out_data.xmat = xmat;

% Final iteration
out_data.istar = i;

% Weight data
out_data.c_mat = c_mat;
out_data.w_mat = w_mat;

% Condition number data
out_data.cond_A_vec = cond_A_vec;


% *************************************************************************
%
% CONTROL SIGNAL
% 
% *************************************************************************

% Initialize empty matrix
umat = zeros(size(tvec,1), 1);

% Calculate control
for k = 1:size(tvec,1)
    
    % Extract state
    x = xmat(k,:)';

    % Extract system variables
    if do_uncertainty
        xs = x(p+1:p+n);        % Known x-dynamics
    else
        xs = x;                 % Known x-dynamics
    end
    
    % Check if in learning phase and set flag appropriately
    if k <= len_tvec_learn
        is_learning = 1;
    else
        is_learning = 0;
    end    
    
    % Evaluate control 
    u = u_t_alg(xs, tvec(k));
                
    % Store control
    umat(k) = u;

end

% Store control signal
out_data.umat = umat;



% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% END MAIN
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CALCULATE DYNAMICS FOR LEARNING PORTION OF RADP ALGORITHM
%
% I.e., steps 1-3 of Algorithm 1 (phase-one and phase-two learning)
%
% NOTE: State partitioning for simulation
%
% In order to perform the RADP algorithm with unmatched uncertainty, the
% following two types of dynamic variables will be necessary for
% simulation:
%
%   System states (including w, x, cf. eqns. (8)-(9))
%   Integrals associated with the least squares minimization (15)
% 
% *************************************************************************
%
% DYNAMIC VARIABLES ASSOCIATED WITH LEAST SQUARES MINIMIZATION (15):
%
% The least squares minimization (15) can be written in the form:
%
%   [A1 A2] [ c_i    = - v
%             w_i ]
%
% Where,
%
% A1 \in R^{l x N1},
%
%       [ \phi_1(x(t_1)) - \phi_1(x(t_0))  ... \phi_{N1}(x(t_1)) - \phi_{N1}(x(t_0))  
% A1 =  [                                  ...
%       [ \phi_1(x(t_l)) - \phi_1(x(t_{l-1}))  ... \phi_{N1}(x(t_l)) - \phi_{N1}(x(t_{l-1}))
% 
% A2 \in R^{l x N2},
%
% \hat{v}_i(t) = u(t) + \Delta(w(t), x(t)) - \hat{u}_i(x(t));
%
%           [ \int_{t_0}^{t_1} \phi_1(x) \hat{v}_i dt  ... \int_{t_0}^{t_1} \phi_{N2}(x) \hat{v}_i dt  
% A2 = 2r * [                                  ...
%           [ \int_{t_{l-1}}^{t_l} \phi_1(x) \hat{v}_i dt  ... \int_{t_{l-1}}^{t_l} \phi_{N2}(x) \hat{v}_i dt
% 
% v \in R^{l}  
%
%       [ \int_{t_0}^{t_1} Q(x) + r \hat{u}_i^2(x) dt  
% v =   [                                  ...
%       [ \int_{t_{l-1}}^{t_l} Q(x) + r \hat{u}_i^2(x) dt
%
%
% Define,
%
% \Phi_{j} : R^n -> R^{Nj}, \Phi_{j}(x) = (\phi_1(x), ..., \phi_{Nj}(x)),
% j = 1, 2, 4
%
% IPsiPsi_mat \in R^{l x N2}
%
%                   [ ( \int_{t_0}^{t_1}  kron(\Phi_{2}(x), \Phi_{2}(x)) )^T    ]
% IPsiPsi_mat =   [                   ...                                       ]
%                   [ ( \int_{t_{l-1}}^{t_l} kron(\Phi_{2}(x), \Phi_{2}(x)) )^T ]
%
% Where kron(A,B) is the Kronecker product of two matrices/vectors A, B
% (see Mathworks for description).
%
% Then, it can be verified algebraically that the following identities
% hold:
%           [ (\int_{t_0}^{t_1} \Phi_{2}(x)(z + \Delta) dt)^T      ]
% A2 = 2r * [                  ...                                  ]
%           [ (\int_{t_{l-1}}^{t_l} \Phi_{2}(x)(z + \Delta) dt)^T  ]
%
%    - 2r * IPsiPsi_mat * kron(w_i, eye(N2)) 
%
%   Where w_i is the current actor NN weight vector, i.e., \hat{u}_i(x) =
%   \sum_{j = 1}^{N2} w_i(j) * \phi_j(x),
% 
%
%       [ \int_{t_0}^{t_1} Q(x) dt      ]
% v =   [                  ...          ]  - r * IPsiPsi_mat * kron(w_i, w_i)
%       [ (\int_{t_{l-1}}^{t_l} Q(x) dt ]
%
%
% The above two identities motivate the need for the dynamic variables:
%
%   IPsiPsi       (N2^2-dim.) Associated with the integral of     
%                   kron(\Phi_{2}(x), \Phi_{2}(x))
%   IPsiupd        (N2-dim.) Associated with the integral of     
%                   \Phi_{2}(x)(u + \Delta)
%   IQ              (1-dim.) Associated with the integral of Q(x)
%
%
% *************************************************************************
%
% FINAL STATE PARTITION:
%
% WITH UNCERTAINTY:
%
% x_sim =   [   w           R^p
%               x           R^n
%               IPsiPsi   R^{N_2^2}
%               IPsiupd    R^{N_2}
%               IQ          R^1       ]
%
% WITHOUT UNCERTAINTY:
%
% x_sim =   [   x           R^n
%               IPsiPsi   R^{N_2^2}
%               IPsiupd    R^{N_2}
%               IQ          R^1       ]
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

function xdot = odefunct(t, x)

% Global variables
global sys;
global Q;
global basis;
global do_uncertainty;

% Get system order
if do_uncertainty
    p = sys.p;          % Order of unknown w-dynamics
end
n = sys.n;          % Order of known x-dynamics

% Extract system variables
if do_uncertainty
    w = x(1:p);                 % Unknown w-dynamics
    xs = x(p+1:p+n);            % Known x-dynamics
else
    xs = x(1:n);                % Known x-dynamics
end


% Total system dynamics
if do_uncertainty
    x_sys = [   w
                xs  ];
else
    x_sys = xs;
end


% ***********************
%       
% EVALUATE SYSTEM DYNAMICS
%
% NOTE: Output data "fx" is partitioned as follows (cf. eval_f.m):
%
% WITH UNCERTAINTY:
%
% fx =  [   Delta_w(w, x)   : w-dynamic uncertainty (R^p)
%           f(x)            : Nominal x-drift dynamics (R^n)
%           Delta(w, x)     : x-dynamic uncertainty (R^1)       ]
%
% WITHOUT UNCERTAINTY:
%
% fx =  [   f(x)            : Nominal x-drift dynamics (R^n)    ]
%

fx = eval_f(x_sys, sys);

% Unpack drift dynamics terms
if do_uncertainty
    Delta_w = fx(1:p);
    f_x = fx(p+1:p+n);
    Delta = fx(p+n+1);
else
    f_x = fx;
    Delta = 0;
end

% Evaluate input gain matrix
gx = eval_g(xs, sys);

% Evaluate basis functions
% Phix = eval_phi(xs, basis.Phi.tag);
Psix = eval_phi(xs, basis.Psi);

% Evaluate Q(x)
Qx = eval_Q(xs, Q);

% Evaluate control signal u(x(t)) = \hat{u}_0(x(t)) + e(t)
u = u_t_alg(xs, t);

% Calculate state derivative \dot{w} (cf. eqn. (8))
if do_uncertainty
    dw = Delta_w;
end

% Calculate state derivative \dot{x} (cf. eqn. (9)) 
dxs = f_x + gx * (u + Delta);

% ***********************
%       
% EVALUATE DYNAMIC VARIABLES ASSOCIATED WITH LEAST SQUARES MINIMIZATION
% (15) (SEE ABOVE)
%

dIPsiPsi = kron(Psix, Psix);
dIPsiupd = Psix * (u + Delta);
dIQ = Qx;


% ***********************
%       
% PACKAGE STATE DERIVATIVE OUTPUT
%
% x_sim =   [   w           R^p
%               x           R^n
%               IPsiPsi   R^{N_2^2}
%               IPsiupd    R^{N_2}
%               IQ          R^1         ]
%
if do_uncertainty
    xdot = [    dw
                dxs
                dIPsiPsi   
                dIPsiupd   
                dIQ             ];
else
    xdot = [    dxs
                dIPsiPsi   
                dIPsiupd   
                dIQ             ];
end


        
%%        
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CALCULATE DYNAMICS FOR LEARNING PORTION OF RADP ALGORITHM
%
% I.e., steps 4-5 of Algorithm 1 (approximate robust optimal control policy
% (31))
%
% NOTE: State partitioning for simulation
%
% The least squares minimization (15) no longer needs to be performed. As
% such, the associated integration variables no longer need to be
% simulated. This portion of the algorithm only requires the system states
% as dynamic variables.
%
% FINAL STATE PARTITION:
%
% WITH UNCERTAINTY:
%
% x_sim =   [   w           R^p
%               x           R^n     ]
%
% WITHOUT UNCERTAINTY:
%
% x_sim =   [   x           R^n     ]
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

function xdot = odefunct_final(t, x)

% Global variables
global sys;
global do_uncertainty;

% Get system order
if do_uncertainty
    p = sys.p;          % Order of unknown w-dynamics
end
n = sys.n;          % Order of known x-dynamics

% Extract system variables
if do_uncertainty
    w = x(1:p);                 % Unknown w-dynamics
    xs = x(p+1:p+n);            % Known x-dynamics
else
    xs = x(1:n);                % Known x-dynamics
end

% Total system dynamics
if do_uncertainty
    x_sys = [   w
                xs  ];
else
    x_sys = xs;
end


% ***********************
%       
% EVALUATE SYSTEM DYNAMICS
%
% NOTE: Output data "fx" is partitioned as follows (cf. eval_f.m):
%
% WITH UNCERTAINTY:
%
% fx =  [   Delta_w(w, x)   : w-dynamic uncertainty (R^p)
%           f(x)            : Nominal x-drift dynamics (R^n)
%           Delta(w, x)     : x-dynamic uncertainty (R^1)       ]
%
% WITHOUT UNCERTAINTY:
%
% fx =  [   f(x)            : Nominal x-drift dynamics (R^n)    ]
%

fx = eval_f(x_sys, sys);

% Unpack drift dynamics terms
if do_uncertainty
    Delta_w = fx(1:p);
    f_x = fx(p+1:p+n);
    Delta = fx(p+n+1);
else
    f_x = fx;
    Delta = 0;
end

% Evaluate input gain matrix
gx = eval_g(xs, sys);

% Evaluate approximate robust optimal control policy u_{ro} (31)
u_ro = u_t_alg(x, t);

% Calculate state derivative \dot{w} (cf. eqn. (8))
if do_uncertainty
    dw = Delta_w;
end

% Calculate state derivative \dot{x} (cf. eqn. (9)) 
dxs = f_x + gx * (u_ro + Delta);
      
        
% ***********************
%       
% PACKAGE STATE DERIVATIVE OUTPUT
%
% x_sim =   [   w           R^p
%               x           R^n    ]
%
if do_uncertainty
    xdot = [    dw
                dxs     ];
else
    xdot = [    dxs     ];
end

        
%%        
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% EVALUATE CONTROL SIGNAL u(t)
%
% The control applied depends on the current stage of the algorithm.
%
% *************************************************************************
% *************************************************************************
% *************************************************************************        

function u = u_t_alg(x, t)

% Global variables
global R;
global basis;
global noise;
global w_i;
global is_learning;
% global tf;
global u_0;
global rho;

global do_uncertainty;

if is_learning
% if t < tf
    
    % *********************************************************************
    %
    % LEARNING PHASE -- INITIAL STABILIZING POLICY u_0(x) (ASSUMPTION 3.2)
    %
	% *********************************************************************
    
    % Evaluate noise
    et = eval_noise(t, 1, noise);
   
    % Evaluate control
    u = eval_u(x, u_0) + et;
    
else

    % *********************************************************************
    %
    % POST-LEARNING -- APPROXIMATE ROBUST OPTIMAL CONTROL POLICY (31)
    %
	% *********************************************************************    
    
    % Evaluate basis functions
    Psix = eval_phi(x, basis.Psi);

    % Evaluate current policy \hat{u}_{i^*+1}
    hat_u_i = w_i' * Psix;

    % If dynamic uncertainty is present, robustify policy
    if do_uncertainty
        
        % Evaluate \rho(||x||^2)
        rhox = eval_rho_radp(norm(x)^2, rho);

        % Evaluate robust redesign policy u_{ro}(x) (cf. eqn. (31))
        % NOTE: The policy u_ro as descirbed in (31) is 
        %
        %   u_ro = (1 + r / 2 * rhox^2) * hat_u_i
        %
        % But as Yu Jiang implemented it in his code, it is
        %
        %   u_ro = (1 + r / 2 * rhox^2) * hat_u_i     +   rhox * hat_u_i
        %

        % As in (31)
        u_ro = (1 + R / 2 * rhox^2) * hat_u_i;       

    %     % As implemented in code
    %     u_ro = (1 + R / 2 * rhox^2 + rhox) * hat_u_i;   
    
        % Control applied is u_{ro}(x)
        u = u_ro;
        
    else
       
        % No dynamic uncertainty. Use PI policy
        u = hat_u_i;
        
    end

end


