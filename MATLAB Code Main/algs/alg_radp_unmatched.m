function out_data = alg_radp_unmatched(alg_settings)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% ROBUST ADAPTIVE DYNAMIC PROGRAMMING (RADP) ALGORITHM -- UNMATCHED
% UNCERTAINTY
%
% Brent Wallace  
%
% 2021-11-06
%
% This program implements the RADP algorithm (with unmatched uncertainty)
% presented in,
%
%   Y. Jiang and Z.-P. Jiang. "Robust adaptive dynamic programming and
%   feedback stabilization of nonlinear systems." IEEE Transactions on
%   Neural Networks and Learning Systems, 25:882-893, 2014.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% out_data = alg_radp_unmatched(alg_settings)
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
%       .Phi1               (Struct) contains critic NN activation function
%                           parameters (cf. eqn. (13)). Has fields:
%           .tag            (String) tag of desired activation functions.
%           .N              (Integer) number of critic NN activation
%                           functions N1.
%       .Phi2               (Struct) contains actor NN activation function
%                           parameters (cf. eqn. (14)). Has fields:
%           .tag            (String) tag of desired activation functions.
%           .N              (Integer) number of actor NN activation
%                           functions N2.
%       .Phi3               (Struct) contains NN activation function
%                           parameters for approximation of \hat{f}_1(x,z)
%                           (cf. eqn. (43)). Has fields:
%           .tag            (String) tag of desired activation functions.
%           .N              (Integer) number of activation functions N3 for
%                           approximation of \hat{f}_1(x,z).
%       .Phi4               (Struct) contains NN activation function
%                           parameters for approximation of \hat{g}_1(x)
%                           (cf. eqn. (44)). Has fields:
%           .tag            (String) tag of desired activation functions.
%           .N              (Integer) number of activation functions N3 for
%                           approximation of \hat{f}_1(x,z).
%   noise                   (Struct) contains info for probing noise. Has
%                           the following fields:
%       .tag                (String) tag of specific probing noise signal
%                           to be injected (see eval_noise.m for options).
%   rho                     (String) tag corresponding to the robust
%                           redesign function \rho(s) (cf. under eqn.
%                           (31)).
%   l                       (Integer) number of samples to collect for each
%                           least-squares minimization (cf. eqns. (47),
%                           (48)).
%   eps1                    (Double) Tolerance \epsilon_1 to determine
%                           phase-1 learning termination (cf. eqn. (59)).
%   eps                     (Double) \epsilon > 0 such that Q(x) -
%                           \epsilon^2 ||x||^2 is a positive-definite
%                           function (cf. Assumption 3.3, Sec. III. B.)
%   tf                      (Double) Length of learning window [0, t_f]
%                           (sec).
%   tsim                    (Double) Length of simulation to run after
%                           learning window (sec). I.e., post-learning
%                           simulation happens [t_f, t_f + tsim].
%   u_0                     (String) Tag corresponding to the initial
%                           stabilizing policy u_0(x,z) (cf. Assumption
%                           4.2). Evaluation of this policy is handled by
%                           eval_u.m, so the policy needs to be encoded
%                           there.
%   w0                      (p-dimensional vector) ICs for unknown dynamics
%                           w(t_0).
%   x0                      (n-dimensional vector) ICs for known x-dynamics
%                           x(t_0).
%   z0                      (Double) ICs for known z-dynamics z(t_0).
%   w_0                     (N2-dimensional vector) ICs for actor NN
%                           weights.
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
%                               x(tvec(k))   R^n
%                               z(tvec(k))   R^1   ]^T
%       .umat               ('simlength' x 1 matrix) Matrix (vector,
%                           really) whose row indexes the time instants
%                           specified in .tvec, and whose 1 column is the
%                           control signal u(t) at the respective time
%                           instant.
%       .istar              (Integer) final termination index of the PI
%                           algorithm.
%       .c_i_mat            ('istar'+1 x N1 matrix) Matrix whose row
%                           indexes the PI index, and whose N1-columns are
%                           the critic NN weights c_i at the respective PI
%                           index i.
%       .w_i_mat            ('istar'+1 x N2 matrix) Matrix whose row
%                           indexes the PI index, and whose N2-columns are
%                           the actor NN weights w_i at the respective PI
%                           index i.
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
global l;
global T;

% Keeps track of if the algorithm is in the learning phase or the final
% phase
global is_learning;

global eps;

global rho;

% % PI index variable
% global i;

% Critic NN parameters
global c_i;

% Actor NN parameters
global w_i;

% Parameters for estimation of \hat{f}_1(x,z)
global w_f;

% Parameters for estimation of \hat{g}_1(x) (cf. eqn. (43))
global w_g;



% *************************************************************************
% 
% UNPACK ALGORITHM SETTINGS/PARAMETERS
% 
% *************************************************************************


preset = alg_settings.preset;               
sys = alg_settings.sys;         % System


p = alg_settings.sys.p;         % Order of unknown w-dynamics
n = alg_settings.sys.n;         % Order of known x-dynamics

Q = alg_settings.Q;
R = alg_settings.R;

noise = alg_settings.noise;

u_0 = alg_settings.u_0;

rho = alg_settings.rho;

l = alg_settings.l;
eps1 = alg_settings.eps1;
eps = alg_settings.eps;

% Simulation times
tf = alg_settings.tf;
tsim = alg_settings.tsim;

% ***********************
%       
% BASIS SETTINGS
%
% See description for details.
%

% Basis struct
basis = alg_settings.basis;

% Critic basis for approximating cost function V_i(x) (cf. eqn. (13))
N1 = alg_settings.basis.Phi1.N;    
% Actor basis for approximating policy u_{i + 1}(x) (cf. eqn. (14))
N2 = alg_settings.basis.Phi2.N;
% Basis for approximating \bar{f_1}(x,z) (cf. eqn. (43))
N3 = alg_settings.basis.Phi3.N;
% Basis for approximating \bar{g_1}(x) (cf. eqn. (44))
N4 = alg_settings.basis.Phi4.N;



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
%   System states (including w, x, z, cf. eqns. (38)-(40))
%   Integrals associated with the least squares minimization (47)
%   Integrals associated with the least squares minimization (48)
%

% Initial conditions
w0 = alg_settings.w0;
xs0 = alg_settings.x0;
z0 = alg_settings.z0; 

% Initial condition for simulation
% See odefunct() for a description of the state partition
%
% x_sim =   [   w           R^p
%               x           R^n
%               z           R^1
%               Iphi2phi2   R^{N_2^2}
%               Iphi2zpd    R^{N_2}
%               IQ          R^1
%               Iphi3       R^{N_3}
%               Iphi4d      R^{N_4}
%               Iupd1       R^1         ]
%
x0_sim = [  w0
            xs0
            z0
            zeros(N2^2, 1)
            zeros(N2, 1);
            0
            zeros(N3, 1)
            zeros(N4, 1)
            0               ];



% ***********************
%       
% MISCELLANEOUS VARIABLES
%

% Initialize critic NN weight vector
c_i = zeros(N1, 1);     % Holds critic NN weights of current iteration

% Initialize actor NN weight vector
w_i = alg_settings.w_0; % Holds actor NN weights of current iteration   

% % Initialize vector to store previous iteration critic NN weights
% c_im1 = zeros(N1, 1);   % Holds critic NN weights of previous iteration
% c_im1(1) = Inf;         % So PI algorithm starts properly

% Initialize vector to store previous iteration actor NN weights
w_im1 = zeros(N2, 1);   % Holds critic NN weights of previous iteration
w_im1(1) = Inf;         % So PI algorithm starts properly

% Initialize weights w_f, w_g (cf. eqns. (43)-(44))
w_f = zeros(N3, 1);
w_g = zeros(N4, 1);

% Sample collection window T
T = tf / l;

% Set learning flag
is_learning = 1;

% Initialize PI index
i = 0;


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
%   Iphi2phi2(k)   R^{N_2^2}
%   Iphi2zpd(k)    R^{N_2}
%   IQ(k)          R^1
%   Iphi3(k)       R^{N_3}
%   Iphi4d(k)      R^{N_4}
%   Iupd1(k)       R^1        
%
%       0 <= k <= l-1
%
% (see odefunct() for a description of these variables) need to be stored
% in corresponding l-row matrices (or l-dimensional vectors) for subsequent
% least squares minimization at the end of the current PI iteration. These
% storage variables are initialized as zeros here.
%

Iphi2phi2_mat = zeros(l,N2^2);
Iphi2zpd_mat = zeros(l,N2);
IQ_vec = zeros(l,1);
Iphi3_mat = zeros(l,N3);
Iphi4d_mat = zeros(l,N4);
Iupd1_vec = zeros(l,1);


% In addition, the non-dynamic variables
%
%   d_Phi1(k) = \Phi_{1}(x(t_{k+1})) - \Phi_{1}(x(t_k))       R^{N_1}
%   d_zeta(k) = \zeta(x(t_{k+1})) - \zeta(x(t_k))               R^1
%
%       0 <= k <= l-1
%
% (see the description of odefunct() for definition of the function
% \Phi_{1}) need to be stored.
%

d_Phi1_mat = zeros(l, N1);
d_zeta_vec = zeros(l, 1);

% An intermediate storage variable used in calculation of d_zeta_vec
d_zeta_0 = zeros(l, 1 + N2);


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
c_i_mat = c_i';

% Stores actor weights w_i at each iteration of the PI algorithm (iteration
% indexed by row)
w_i_mat = w_i';

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
% See steps 1-4 of Algorithm 2, Jiang and Jiang
% 
% *************************************************************************
% *************************************************************************


% *************************************************************************
%
% COLLECT SIMULATION DATA
%
% This loop collects the l data points needed to perform the least-squares
% minimizations (47) and (48).
%  
% *************************************************************************

for k = 0:l-1

    % ***********************
    %       
    % RUN SIMULATION
    %

    tspan = k * T + [0, T];

    [t, x] = ode45(@odefunct, tspan, x0_sim);

    % DEBUGGING: Plot state trajectory
    figure(100)
    hold on
    plot(t,x(:,1:p+n+1));
    title('State Trajectory')
    grid on
    xlabel('Time (sec)')
    ylabel('x(t)')

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
    %
    % x_sim =   [   w           R^p
    %               x           R^n
    %               z           R^1
    %               Iphi2phi2   R^{N_2^2}
    %               Iphi2zpd    R^{N_2}
    %               IQ          R^1
    %               Iphi3       R^{N_3}
    %               Iphi4d      R^{N_4}
    %               Iupd1       R^1         ]

    % State dynamic variables
    w0 = x0(1:p);               % Unknown w-dynamics
    xs0 = x0(p+1:p+n);        % Known x-dynamics
    z0 = x0(p+n+1);           % Known z-dynamics

    % Unpack state data at end of simulation. See odefunct() for a
    % description of the state partition.
    %
    % x_sim =   [   w           R^p
    %               x           R^n
    %               z           R^1
    %               Iphi2phi2   R^{N_2^2}
    %               Iphi2zpd    R^{N_2}
    %               IQ          R^1
    %               Iphi3       R^{N_3}
    %               Iphi4d      R^{N_4}
    %               Iupd1       R^1         ]

    % State dynamic variables
    w1 = x1(1:p);               % Unknown w-dynamics
    xs1 = x1(p+1:p+n);        % Known x-dynamics
    z1 = x1(p+n+1);           % Known z-dynamics

    % Dynamic variables associated with least squares minimizations
    x1_lsq = x1(p+n+2:end);     

    Iphi2phi2 = x1_lsq(1:N2^2);
    Iphi2zpd = x1_lsq(N2^2+1:N2^2+N2);
    IQ = x1_lsq(N2^2+N2+1);
    Iphi3 = x1_lsq(N2^2+N2+1+1:N2^2+N2+1+N3);
    Iphi4d = x1_lsq(N2^2+N2+1+N3+1:N2^2+N2+1+N3+N4);
    Iupd1 = x1_lsq(N2^2+N2+1+N3+N4+1);

    % ***********************
    %       
    % DATA STORAGE -- SYSTEM STATE DATA, DYNAMIC LEAST-SQUARES DATA
    % 

    % Store time data
    tvec = [    tvec
                t       ];

    % Store system state data
    xmat =     [   xmat
                    x(:,1:p+n+1)    ];

    % Store dynamic variables associated with least squares data 

    Iphi2phi2_mat(k+1,:) = Iphi2phi2';
    Iphi2zpd_mat(k+1,:) = Iphi2zpd';
    IQ_vec(k+1) = IQ;

    Iphi3_mat(k+1,:) = Iphi3';
    Iphi4d_mat(k+1,:) = Iphi4d';
    Iupd1_vec(k+1) = Iupd1;

    % ***********************
    %       
    % DATA STORAGE -- REMAINING LEAST SQUARES DATA
    %         

    % Evaluate critic basis functions \Phi_{1}(x) = (\phi_1(x), ...,
    % \phi_{N1}(x)) at x = x(t_k) and x = x(t_{k+1})
    Phi1x_tk = eval_phi(xs0, basis.Phi1);
    Phi1x_tkp1 = eval_phi(xs1, basis.Phi1); 

    % Calculate and store \Phi_{1}(x(t_{k+1})) - \Phi_{1}(x(t_k))
    d_Phi1_mat(k+1,:) = (Phi1x_tkp1 - Phi1x_tk)';

    % Evaluate critic basis functions \Phi_{2}(x) = (\phi_2(x), ...,
    % \phi_{N2}(x)) at x = x(t_k) and x = x(t_{k+1})
    Phi2x_tk = eval_phi(xs0, basis.Phi2);
    Phi2x_tkp1 = eval_phi(xs1, basis.Phi2);     
    
%     % Evaluate policy       
%     hat_u_i_tk = w_i' * Phi2x_tk;
%     hat_u_i_tkp1 = w_i' * Phi2x_tkp1;

    % Evaluate \rho(||x||^2)
    rhox0 = eval_rho_radp(norm(xs0)^2, preset);
    rhox1 = eval_rho_radp(norm(xs1)^2, preset);

    % Evaluate intermediate term d_zeta_0 for later calculation of
    % d_zeta_vec
    %
    % d_zeta_0 \in R^{l x (1 + N2)},
    %
    %            [ z(t_1) - z(t_0)      | (1 + r / 2 * \rho^2(||x(t_1)||^2)) *
    %            [                      |  \Phi_{2}(x(t_1))^T
    %            [                      | -
    %            [                      | (1 + r / 2 * \rho^2(||x(t_0)||^2)) *
    %            [                      | \Phi_{2}(x(t_0))^T
    % d_zeta_0 = [  ...                 |       ...
    %            [ z(t_l) - z(t_{l-1})  | (1 + r / 2 * \rho^2(||x(t_l)||^2)) *
    %            [                      |  \Phi_{2}(x(t_l))^T
    %            [                      | -
    %            [                      | (1 + r / 2 * \rho^2(||x(t_{l-1})||^2)) *
    %            [                      | \Phi_{2}(x(t_{l-1}))^T              
    %
    % NOTE: After PI algorithm terminates, we have
    %   
    %   d_zeta_mat = d_zeta_0 * [ 1 
    %                             -w_i ] 
    %
    
    d_zeta_0(k+1,1) = z1 - z0;
    d_zeta_0(k+1,2:end) = (1 + R / 2 * rhox1^2) * Phi2x_tkp1'...
                            - (1 + R / 2 * rhox0^2) * Phi2x_tk';
    
%     % Evaluate robust redesign policy u_{ro}(x) (cf. eqn. (31))
%     u_ro_tk = (1 + R / 2 * rhox0^2) * hat_u_i_tk;
%     u_ro_tkp1 = (1 + R / 2 * rhox1^2) * hat_u_i_tkp1;
% 
%     % Evaluate state transformation \zeta = z - u_ro (cf. beginning of
%     % Sec. IV. A.)
%     zeta_tk = z0 - u_ro_tk;
%     zeta_tkp1 = z1 - u_ro_tkp1;
% 
%     % Calculate and store \zeta(t_{k+1}) - \zeta(t_k)
%     d_zeta_vec(k+1) = zeta_tkp1 - zeta_tk;

    % ***********************
    %       
    % PREPARE FOR NEXT SIMULATION
    % 

    % IC for next simulation. System state variables for next
    % simulation are carried over from their final values in this
    % simulation. Integration variables associated with least squares
    % minimizations are reset to zero.  
    %
    % See odefunct() for a description of the state partition
    %
    % x_sim =   [   w           R^p
    %               x           R^n
    %               z           R^1
    %               Iphi2phi2   R^{N_2^2}
    %               Iphi2zpd    R^{N_2}
    %               IQ          R^1
    %               Iphi3       R^{N_3}
    %               Iphi4d      R^{N_4}
    %               Iupd1       R^1         ]
    %
    x0_sim = [  w1
                xs1
                z1
                zeros(N2^2, 1)
                zeros(N2, 1);
                0
                zeros(N3, 1)
                zeros(N4, 1)
                0               ];      


end


% *************************************************************************
%
% PHASE 1 LEARNING: PI LOOP
%
% NOTE: As per eqn. (59), the termination condition of the PI algorithm is
%
%   || \hat{c}_i - \hat{c}_{i-1} ||^2 <= \epsilon_1
%
% I.e., termination depends on convergence of the critic NN weights.
% But as Yu Jiang has implemented in his code, the termination condition is
%
%   || \hat{w}_i - \hat{w}_{i-1} || <= \epsilon_1
%
% I.e., termination depends on convergence of the actor NN weights, and on
% the FIRST power of the norm.
%
% We have followed his code implementation for termination here (i.e.,
% actor weights are checked).
%
% *************************************************************************

while norm(w_i - w_im1) > eps1
    
    
    % *********************************************************************
    %       
    % LEAST SQUARES MINIMIZATION (47)
    %
    % Defined as in the description of odefunct(), the least squares
    % minimization (47) has the form:
    %
    %   [A1 A2] [ c_i    = - v
    %             w_i ]
    %
    % Expressed in terms of the variables we have declared above, this is:
    %
    % A1 = d_Phi1_mat;
    % A2 = 2r * ( Iphi2zpd_mat - Iphi2phi2_mat * kron(w_i, eye(N2)) )
    % v = IQ_vec + Iphi2phi2_mat * kron(w_i, w_i)
    %
    
    % Least squares variables
    A1 = d_Phi1_mat;
    A2 = 2 * R * ( Iphi2zpd_mat - Iphi2phi2_mat * kron(w_i, eye(N2)) );
    v = IQ_vec + Iphi2phi2_mat * kron(w_i, w_i);
    
    % Perform least squares
    lsq = [ A1, A2 ]\(-v);
    
    % Store current w_i as w_{i-1}, extract new least squares solutions
%     c_im1 = c_i;
    w_im1 = w_i;
    c_i = lsq(1:N1);
    w_i = lsq(N1+1:end);
    
    % Store c_i, w_i for later plotting
    c_i_mat =   [   c_i_mat
                    c_i'        ];
    w_i_mat =   [   w_i_mat
                    w_i'        ];
                
    % *********************************************************************
    %       
    % PREPARE FOR NEXT PI ITERATION
    %    
    
    % Increment PI index
    i = i + 1;
    
    
end


% *************************************************************************
%
% PHASE 2 LEARNING: PERFORM LEAST SQUARES MINIMIZATION (48)
%
% Defined as in the description of odefunct(), the least squares
% minimization (48) (i.e., as implemented in Yu Jiang's code) has the form:
%
%   [B1 B2] [ w_f    = Z - U
%             w_g ]
%
% Expressed in terms of the variables we have declared above, this is:
%
% B1 = Iphi3_mat
% B2 = Iphi4d_mat
% Z = d_zeta_vec
% U = Iupd1_vec
%
% *************************************************************************

% Calculation of d_zeta_vec
d_zeta_vec = d_zeta_0 * [ 1 ; -w_i];

% Least squares variables
B1 = Iphi3_mat;
B2 = Iphi4d_mat;
Z = d_zeta_vec;
U = Iupd1_vec;

% Perform least squares
lsq = [ B1, B2 ]\(Z - U);

% Extract least-squares optimal weights w_f, w_g
w_f = lsq(1:N3);
w_g = lsq(N3+1:end);

% DEBUGGING: Print final weights
c_i
w_i
w_f
w_g


%%
% *************************************************************************
% *************************************************************************
%
% APPLY APPROXIMATE ROBUST OPTIMAL CONTROL POLICY (54)
%
% See steps 5-6 of Algorithm 2, Jiang and Jiang
% 
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
% x_sim =   [   w           R^p
%               x           R^n
%               z           R^1      ]
%

x0_sim =    [   w1
                xs1
                z1      ];
            
% Time span for simulation
tspan = [tf, tf + tsim];

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
xmat =     [   xmat
               x(:,1:p+n+1)    ];

            
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


out_data.tvec = tvec;
out_data.xmat = xmat;

out_data.istar = i;
out_data.c_i_mat = c_i_mat;
out_data.w_i_mat = w_i_mat;




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
%     w = x(1:p);             % Unknown w-dynamics
    xs = x(p+1:p+n);        % Known x-dynamics
    z = x(p+n+1);           % Known z-dynamics

    % Known system dynamics
    x_known =   [   xs
                    z   ];
                
    % Check if in learning phase and set flag appropriately
    if k <= len_tvec_learn
        is_learning = 1;
    else
        is_learning = 0;
    end
    
    % Evaluate control 
    u = u_t_alg(x_known, tvec(k));
                
    % Store control
    umat(k) = u;

end

% *************************************************************************
%
% STORE CONTROL SIGNAL
% 

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
% I.e., steps 1-4 of Algorithm 2 (phase-one and phase-two learning)
%
% NOTE: State partitioning for simulation
%
% In order to perform the RADP algorithm with unmatched uncertainty, the
% following three types of dynamic variables will be necessary for
% simulation:
%
%   System states (including w, x, z, cf. eqns. (38)-(40))
%   Integrals associated with the least squares minimization (47)
%   Integrals associated with the least squares minimization (48)
% 
% *************************************************************************
%
% DYNAMIC VARIABLES ASSOCIATED WITH LEAST SQUARES MINIMIZATION (47):
%
% The least squares minimization (47) can be written in the form:
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
% \tilde{v}_i(t) = z(t) + \Delta(w(t), x(t)) - \hat{u}_i(x(t));
%
%           [ \int_{t_0}^{t_1} \phi_1(x) \tilde{v}_i dt  ... \int_{t_0}^{t_1} \phi_{N2}(x) \tilde{v}_i dt  
% A2 = 2r * [                                  ...
%           [ \int_{t_{l-1}}^{t_l} \phi_1(x) \tilde{v}_i dt  ... \int_{t_{l-1}}^{t_l} \phi_{N2}(x) \tilde{v}_i dt
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
% Iphi2phi2_mat \in R^{l x N2}
%
%                   [ ( \int_{t_0}^{t_1}  kron(\Phi_{2}(x), \Phi_{2}(x)) )^T    ]
% Iphi2phi2_mat =   [                   ...                                       ]
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
%    - 2r * Iphi2phi2_mat * kron(w_i, eye(N2)) 
%
%   Where w_i is the current actor NN weight vector, i.e., \hat{u}_i(x) =
%   \sum_{j = 1}^{N2} w_i(j) * \phi_j(x),
% 
%
%       [ \int_{t_0}^{t_1} Q(x) dt      ]
% v =   [                  ...          ]  - r * Iphi2phi2_mat * kron(w_i, w_i)
%       [ (\int_{t_{l-1}}^{t_l} Q(x) dt ]
%
%
% The above two identities motivate the need for the dynamic variables:
%
%   Iphi2phi2       (N2^2-dim.) Associated with the integral of     
%                   kron(\Phi_{2}(x), \Phi_{2}(x))
%   Iphi2zpd        (N2-dim.) Associated with the integral of     
%                   \Phi_{2}(x)(z + \Delta)
%   IQ              (1-dim.) Associated with the integral of Q(x)
%
% *************************************************************************
%
% DYNAMIC VARIABLES ASSOCIATED WITH LEAST SQUARES MINIMIZATION (48):
%
% There are significant discrepancies between the least-squares
% minimization (48), and how the least-squares minimization (48) is
% implemented in Yu Jiang's code. In Jiang's code, (48) reads:
%
%   \zeta(t_{k+1}) - \zeta(t_{k})
%       = \int_{t_k}^{t_{k+1}} [ \sum_{j=1}^{N3} \hat_{w}_{f,j} \psi_j(x,z)
%                                - \sum_{j=1}^{N4-1} \hat_{w}_{g,j}
%                                \phi_j(x) ] dt
%       + \int_{t_k}^{t_{k+1}} (u + \Delta_1) dt + \tilde{e}_k,
%
%   k = 0, ..., l-1
%
% As one can see, this is quite different from the minimization (48)
% depicted in the paper. Specifically, the difference on the LHS is taken
% as a difference of linear terms in \zeta. Also, a \zeta multuplication
% term is missing in both integrands of the code implementation, whereas it
% is seen in (48) in the paper.
%
% THE CODE HERE IMPLEMENTS THE LEAST-SQUARES MINIMIZATION AS IMPLEMENTED IN
% YU JIANG'S CODE; I.E., AS DEPICTED IN THE ABOVE EQUATION.
%
% The above least-squares may be written,
%
%   [B1 B2] [ w_f    = Z - U
%             w_g ]
%
% Where,
%
% B1 \in R^{l x N3},
%
%      [ \int_{t_0}^{t_1} \psi_1(x,z) dt  ... \int_{t_0}^{t_1} \psi_{N3}(x,z) dt  
% B1 = [                                  ...
%      [ \int_{t_{l-1}}^{t_l} \psi_1(x,z) dt  ... \int_{t_{l-1}}^{t_l} \psi_{N3}(x,z) dt
% 
% B2 \in R^{l x N4},
%
%      [ \int_{t_0}^{t_1} \phi_0(x) \Delta dt  ... \int_{t_0}^{t_1} \phi_{N4-1}(x) \Delta dt  
% B2 = [                                  ...
%      [ \int_{t_{l-1}}^{t_l} \phi_0(x) \Delta dt  ... \int_{t_{l-1}}^{t_l} \phi_{N4-1}(x) \Delta dt
% 
% Z \in R^{l}  
%
%       [ \zeta(t_1) - \zeta(t_0)       ]
% Z =   [           ...                 ]
%       [ \zeta(t_{l}) - \zeta(t_{l-1}) ]
%
% U \in R^{l}  
%
%       [ \int_{t_0}^{t_1} (u + \Delta_1) dt     ]
% U =   [           ...                         ]
%       [ \int_{t_{l-1}}^{t_l} (u + \Delta_1) dt ]
%
% Define,
%
% \Phi_{j} : R^n -> R^{Nj}, \Phi_{j}(x) = (\phi_1(x), ..., \phi_{Nj}(x)),
% j = 1, 2, 4
%
% \psi_{c} : R^n x R -> R^{N3}, 
%   \psi_{c}(x, z) = (\psi_1(x,z), ..., \psi_{N3}(x,z))
%
% The above motivates the need for the following dynamic variables:
%
%   Iphi3           (N3-dim.) Associated with the integral of \psi_{c}(x,z)
%   Iphi4d          (N4-dim.) Associated with the integral of     
%                   \Phi_{4}(x) \Delta
%   Iupd1           (1-dim.) Associated with the integral of u + \Delta_1
%
% *************************************************************************
%
% FINAL STATE PARTITION:
%
% x_sim =   [   w           R^p
%               x           R^n
%               z           R^1
%               Iphi2phi2   R^{N_2^2}
%               Iphi2zpd    R^{N_2}
%               IQ          R^1
%               Iphi3       R^{N_3}
%               Iphi4d      R^{N_4}
%               Iupd1       R^1         ]
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

function xdot = odefunct(t, x)

% Global variables
global sys;
global Q;
global basis;

% Get system order
p = sys.p;          % Order of unknown w-dynamics
n = sys.n;          % Order of known x-dynamics

% Extract system variables
w = x(1:p);             % Unknown w-dynamics
xs = x(p+1:p+n);        % Known x-dynamics
z = x(p+n+1);           % Known z-dynamics

% Known system dynamics
x_known =   [   xs
                z   ];

% Total system dynamics
x_sys = [   w
            xs
            z   ];
        
% ***********************
%       
% EVALUATE SYSTEM DYNAMICS
%
% NOTE: Output data "fx" is partitioned as follows (cf. eval_f.m):
%
% fx =  [   Delta_w(w, x)   : w-dynamic uncertainty (R^p)
%           f(x)            : Nominal x-drift dynamics (R^n)
%           Delta(w, x)     : x-dynamic uncertainty (R^1)
%           f_1(x, z)       : Nominal z-drift dynamics (R^1)
%           Delta_1(w, x, z): z-dynamic uncertainty (R^1)       ]
%

fx = eval_f(x_sys, sys);

% Unpack drift dynamics terms
Delta_w = fx(1:p);
f_x = fx(p+1:p+n);
Delta = fx(p+n+1);
f_1 = fx(p+n+1+1);
Delta_1 = fx(end);

% Evaluate input gain matrix
gx = eval_g(xs, sys);

% Evaluate basis functions
% Phi1x = eval_phi(xs, basis.Phi1.tag);
Phi2x = eval_phi(xs, basis.Phi2);
Phi3xz = eval_phi(x_known, basis.Phi3);
Phi4x = eval_phi(xs, basis.Phi4);

% Evaluate Q(x)
Qx = eval_Q(xs, Q);

% Evaluate control signal u(x(t)) = \hat{u}_0(x(t)) + e(t)
u = u_t_alg(x_known, t);

% Calculate state derivative \dot{w} (cf. eqn. (38))
dw = Delta_w;

% Calculate state derivative \dot{x} (cf. eqn. (39)) 
dxs = f_x + gx * (z + Delta);

% Calculate state derivative \dot{z} (cf. eqn. (40))
dz = f_1 + u + Delta_1;

% ***********************
%       
% EVALUATE DYNAMIC VARIABLES ASSOCIATED WITH LEAST SQUARES MINIMIZATION
% (47) (SEE ABOVE)
%

dIphi2phi2 = kron(Phi2x, Phi2x);
dIphi2zpd = Phi2x * (z + Delta);
dIQ = Qx;

% ***********************
%       
% EVALUATE DYNAMIC VARIABLES ASSOCIATED WITH LEAST SQUARES MINIMIZATION
% (48) (SEE ABOVE)
%

dIphi3 = Phi3xz;
dIphi4d = Phi4x * Delta;
dIupd1 = u + Delta_1;



% ***********************
%       
% PACKAGE STATE DERIVATIVE OUTPUT
%
% x_sim =   [   w           R^p
%               x           R^n
%               z           R^1
%               Iphi2phi2   R^{N_2^2}
%               Iphi2zpd    R^{N_2}
%               IQ          R^1
%               Iphi3       R^{N_3}
%               Iphi4d      R^{N_4}
%               Iupd1       R^1         ]
%
xdot = [    dw
            dxs
            dz
            dIphi2phi2   
            dIphi2zpd   
            dIQ         
            dIphi3       
            dIphi4d     
            dIupd1          ];

        
%%        
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CALCULATE DYNAMICS FOR LEARNING PORTION OF RADP ALGORITHM
%
% I.e., steps 5-6 of Algorithm 2 (approximate robust optimal control policy
% (54))
%
% NOTE: State partitioning for simulation
%
% The least squares minimizations (47), (48) no longer need to be
% performed. As such, the associated integration variables no longer need
% to be simulated. This portion of the algorithm only requires the system
% states as dynamic variables.
%
% FINAL STATE PARTITION:
%
% x_sim =   [   w           R^p
%               x           R^n
%               z           R^1      ]
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

function xdot = odefunct_final(t, x)

% Global variables
global sys;


% Get system order
p = sys.p;          % Order of unknown w-dynamics
n = sys.n;          % Order of known x-dynamics

% Extract system variables
w = x(1:p);             % Unknown w-dynamics
xs = x(p+1:p+n);        % Known x-dynamics
z = x(p+n+1);           % Known z-dynamics

% Known system dynamics
x_known =   [   xs
                z   ];

% Total system dynamics
x_sys = [   w
            xs
            z   ];
        
% ***********************
%       
% EVALUATE SYSTEM DYNAMICS
%
% NOTE: Output data "fx" is partitioned as follows (cf. eval_f.m):
%
% fx =  [   Delta_w(w, x)   : w-dynamic uncertainty (R^p)
%           f(x)            : Nominal x-drift dynamics (R^n)
%           Delta(w, x)     : x-dynamic uncertainty (R^1)
%           f_1(x, z)       : Nominal z-drift dynamics (R^1)
%           Delta_1(w, x, z): z-dynamic uncertainty (R^1)       ]
%

fx = eval_f(x_sys, sys);

% Unpack drift dynamics terms
Delta_w = fx(1:p);
f_x = fx(p+1:p+n);
Delta = fx(p+n+1);
f_1 = fx(p+n+1+1);
Delta_1 = fx(end);

% Evaluate input gain matrix
gx = eval_g(xs, sys);

% Evaluate approximate robust optimal control policy u_{ro1} (54)
u_ro1 = u_t_alg(x_known, t);

% Calculate state derivative \dot{w} (cf. eqn. (38))
dw = Delta_w;

% Calculate state derivative \dot{x} (cf. eqn. (39)) 
dxs = f_x + gx * (z + Delta);

% Calculate state derivative \dot{z} (cf. eqn. (40))
dz = f_1 + u_ro1 + Delta_1;    
    

% ***********************
%       
% PACKAGE STATE DERIVATIVE OUTPUT
%
% x_sim =   [   w           R^p
%               x           R^n
%               z           R^1   ]
%
xdot = [    dw
            dxs
            dz      ];
        
        
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

global w_f;
global w_g;

if is_learning
% if t < tf
    
    % *********************************************************************
    %
    % LEARNING PHASE -- INITIAL STABILIZING POLICY u_0(x,z) (ASSUMPTION
    % 4.2)
    %
	% *********************************************************************
       
    % Evaluate noise
    et = eval_noise(t, 1, noise);
   
    % Evaluate control
    u = eval_u(x, u_0) + et;
    
else

    % *********************************************************************
    %
    % POST-LEARNING -- APPROXIMATE ROBUST OPTIMAL CONTROL POLICY u_{r01}(x)
    % (54)
    %
	% *********************************************************************    
    
    % Extract state data
    % NOTE: Input is x_known = [ xs ; z ]
    xs = x(1:end-1);
    z = x(end);
      
    % Evaluate basis functions
    % Phi1x = eval_phi(xs, basis.Phi1.tag);
    Phi2x = eval_phi(xs, basis.Phi2);
    Phi3xz = eval_phi(x, basis.Phi3);
    Phi4x = eval_phi(xs, basis.Phi4);

    % Evaluate current policy \hat{u}_{i^*+1}
    hat_u_i = w_i' * Phi2x;

    % Evaluate \hat{f}_1(x,z) (cf. eqn. (43))
    hat_f_1 = w_f' * Phi3xz;

    % Evaluate \hat{g}_1(x) (cf. eqn. (44))
    hat_g_1 = w_g' * Phi4x;

    % Evaluate \rho(||x||^2)
    % NOTE: \rho_1(s) = 2 \rho(s/2) (cf. under eqn. (54))
    rhox = eval_rho_radp(norm(xs)^2, rho);
    rho1x = 2 * eval_rho_radp(0.5 * norm(xs)^2, rho);

    % Evaluate robust redesign policy u_{ro}(x) (cf. eqn. (31))
    % NOTE: The policy u_ro as descirbed in (31) is 
    %
    %   u_ro = (1 + r / 2 * rhox^2) * hat_u_i
    %
    % But as Yu Jiang implemented it in his code, it is
    %
    %   u_ro = (1 + r / 2 * rhox^2) * hat_u_i     +   rhox * hat_u_i
    %

%     % As in (31)
%     u_ro = (1 + R / 2 * rhox^2) * hat_u_i;     

    % As implemented in code
    u_ro = (1 + R / 2 * rhox^2 + rhox) * hat_u_i;   

    % Evaluate state transformation \zeta = z - u_ro (cf. beginning of Sec.
    % IV. A.)
    zeta = z - u_ro;

    % Evaluate the term \rho_1(||X_1||^2) (cf. eqn. (54))
    % NOTE: X_1 = [x^T, \zeta]^T (cf. under eqn. (54))
    X_1 =   [   x
                zeta    ];

    rho1_X_1 = 2 * eval_rho_radp(0.5 * norm(X_1)^2, rho);    

    % Evaluate the terms \rho(zeta^2) (cf. eqn. (54)), \rho_1(zeta^2)
    rho_zeta = eval_rho_radp(zeta^2, rho); 
    rho1_zeta = 2 * eval_rho_radp(0.5 * zeta^2, rho);

    % Evaluate the approximate robust optimal control policy u_{ro1} (54)
    %
    % NOTE: The policy u_{ro1} as shown in (54) varies from how it was
    % implemented in Yu Jiang's code. Both versions are coded below.
    %
    % List of differences:
    %   Line 2:     '+' in (54), vs. '-' in code
    %   Line 5:     '1/4' in (54), vs. '1/16' in code
    %   Line 6:     '\rho(.)' in (54) vs. '\rho_1(.) in code ( in both
    %               numerator and denominator terms)
    %

    % % As in (54):
    % u_ro1 = - hat_f_1 ...
    %         + 2 * R * hat_u_i ...
    %         - 0.25 * hat_g_1^2 * rho1_X_1^2 * zeta ...
    %         - eps^2 * zeta ...
    %         - 1 / 4 * rho1_X_1^2 * zeta ...
    %         - 0.5 * eps^2 * rho_zeta^2 * zeta / rhox^2;

    % As implemented in Yu Jiang's code:
    u_ro1 = - hat_f_1 ...
            - 2 * R * hat_u_i ...
            - 0.25 * hat_g_1^2 * rho1_X_1^2 * zeta ...
            - eps^2 * zeta ...
            - 1 / 16 * rho1_X_1^2 * zeta ...
            - 0.5 * eps^2 * rho1_zeta^2 * zeta / rho1x^2;
    
    % Control applied is u_{ro1}(x)
    u = u_ro1;

end

        

