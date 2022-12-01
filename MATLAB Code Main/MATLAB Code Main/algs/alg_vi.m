function out_data = alg_vi(alg_settings)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% VALUE ITERATION (VI) ALGORITHM
%
% Brent Wallace  
%
% 2022-01-06
%
% This program implements the VI algorithm presented in,
%
%   T. Bian and Z.-P. Jiang. "Reinforcement Learning and Adaptive Optimal
%   Control for Continuous-Time Nonlinear Systems: A Value Iteration
%   Approach." IEEE Transactions on Neural Networks and Learning Systems,
%   Accepted for Publication, 2021.
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
%       Phi                 (Struct) Contains parameters for the critic NN
%                           basis. Has the following fields:
%           tag             (String) Tag of desired activation functions.
%           N               (Integer) Number of critic activation functions
%                           N_1.
%       Psi                 (Struct) Contains parameters for the
%                           Hamiltonian NN basis which are associated with
%                           the actor. Has the following fields:
%           tag             (String) Tag of desired activation functions.
%           N               (Integer) Number of activation functions of x
%                           and u N_2. 
%       Theta               (Struct) Contains parameters for the remaining
%                           activations of the Hamiltonian NN basis (i.e.,
%                           those which do not pertain to the actor). Has
%                           the following fields:
%           tag             (String) Tag of desired activation functions.
%           N               (Integer) Number of activation functions of x
%                           only N_3.
%   noise                   (Struct) contains info for probing noise. Has
%                           the following fields:
%       tag                 (String) tag of specific probing noise signal
%                           to be injected (see eval_noise.m for options).
%   u_0                     (String, optional) Tag corresponding to an
%                           initial policy u_0(x). NOTE: Initial policy not
%                           required for VI algorithm. If one is not used,
%                           then simply declare this vield as '0'.
%   sf                      (Double) Amount of time to run weight
%                           differential equation learning for (sec).
%   tf                      (Double) Length of learning window [0, t_f]
%                           (sec).
%   tsim                    (Double) Length of simulation to run after
%                           learning window (sec). I.e., post-learning
%                           simulation happens [t_f, t_f + tsim].
%   maxstep                 (Double) Max step size for ode45 (sec).
%   r1_R0                   (Bool) 1 = last basis function in Hamiltonian
%                           basis is r(x,u) = Q(x) + u^T R u. 0 = last
%                           basis function is u^T R u. NOTE: If r(x,u) is
%                           the last basis function, basis functions for
%                           Q(x) must be included in the \Theta(x) basis.
%   include_in_Sigma        (Bool) 1 = include last basis function
%                           explicitly in the Hamiltonian basis. 0 = don't
%                           include explicitly in the basis.
%   int_H_ode45_1_man_0     (Bool) 1 = perform Hamiltonian function
%                           integral for the critic weight update (13) via
%                           ode45 (=1, slower, more numerically reliable)
%                           or manually with pre-existing state trajectory
%                           data (=0, faster).
%   x0                      (n-dimensional vector) ICs x(t_0).
%   c_0                     (N1-dimensional vector) ICs for critic NN
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
%   tvec                    ('simlength'-dimensional vector) vector of
%                           'simlength' time indices corresponding to the
%                           simulation time instants over the course of the
%                           algorithm execution.
%   xmat                    ('simlength' x (p+n+1) matrix) Matrix whose row
%                           indexes the time instants specified in .tvec,
%                           and whose n-columns are the state vector at the
%                           respective time instant.
%   umat                    ('simlength' x m matrix) Matrix (vector,
%                           really) whose row indexes the time instants
%                           specified in .tvec, and whose m columns are the
%                           control signal u(t) at the respective time
%                           instant.
%   lentvec_learn           (Integer) Number of simulation data in the
%                           learning phase of the algorithm; i.e., over the
%                           interval [0, t_f].
%   svec                    ('max_iter'-dim vector) Vector whose ith entry
%                           is the ith time instant of the VI learning
%                           dynamic weight update. So svec(1) = 0,
%                           svec(max_iter) = s_f.
%   c_mat                   ('max_iter' x N1 matrix) Matrix whose ith row
%                           indexes the ith VI time index svec(i), and
%                           whose N1-columns are the critic NN weights c_s
%                           at the respective index svec(i).
%   w_mat                   ('max_iter' x N2 matrix) Matrix whose ith row
%                           indexes the ith VI time index svec(i), and
%                           whose N2-columns are the Hamiltonian NN weights
%                           w_s (i.e., the weights associated with the
%                           actor) at the respective index svec(i).
%   v_mat                   ('max_iter' x N3 matrix) Matrix whose ith row
%                           indexes the ith VI time index svec(i), and
%                           whose N3-columns are the Hamiltonian NN weights
%                           v_s (i.e, the weights NOT associated with the
%                           actor) at the respective index svec(i).
%   cond_A_vec              (2-dim. vector) The first index is the
%                           condition number cond(K_{\phi}(t_f)). The
%                           second is the condition number
%                           cond(K_{\sigma}(t_f))
%   num_iter                (Integer) Number of iterations performed by the
%                           algorithm before termination.
%     
% *************************************************************************
%
% REMARK ON NOTATION
%
% *************************************************************************
%
% In T. Bian, Z.-P. Jiang, the there is the notation
%
% PI INDEX VARIABLE: k
%
% CRITIC NN:
%
%   ACTIVATION FUNCTIONS:
%   {\phi_i(x))}_{i=1}^{\infty}
%
%   TRUNCATION:
%   {\phi_i(x))}_{i=1}^{N_1}
%
%   WEIGHTS:
%   w_k \in R^{N_1}
%
% HAMILTONIAN NN:
%
%   ACTIVATION FUNCTIONS:
%   {\psi_i(x,u)}_{i=1}^{\infty} = 
%
%       {\psi_i^{(0)}(x)}_{i=1}^{\infty}
%       U
%       {\psi_i^{(1)}(x) * u}_{i=1}^{\infty}
%       U
%       {u^T R u}
%
%   TRUNCATION:
%   {\psi_i^{(0)}(x)}_{i=1}^{N_{2,0}}
%   {\psi_i^{(1)}(x)}_{i=1}^{N_{2,1}}
%   N_2 = N_{2,0} + N_{2,1}
%
%   \Phi(x) = (\phi_1(x), ... , \phi_{N1}(x))
%   \Psi(x,u) = [       \psi_1^{(0)}(x)
%                       ...
%                       \psi_{N_{2,0}}^{(0)}(x)
%                       \psi_1^{(1)}(x) u
%                       ...
%                       \psi_{N_{2,1}}^{(1)}(x) u
%                       u^T R u             ]
%
%   WEIGHTS:
%   c_k^{(0)} \in R^{N_{2,0}}
%   c_k^{(1)} \in R^{N_{2,1}}
%
% Unfortunately, none of this notation is consistent with the standard
% notation we have adopted. We thus make the name changes:
%
% PI INDEX VARIABLE: i
%
% CRITIC NN:
%
%   ACTIVATION FUNCTIONS:
%   {\phi_j(x))}_{j=1}^{\infty}
%
%   TRUNCATION:
%   {\phi_j(x))}_{j=1}^{N_1}
%
%   WEIGHTS:
%   w_s \in R^{N_1}
%
% HAMILTONIAN NN:
%
%   ACTIVATION FUNCTIONS:
%   {\sigma_j(x,u)}_{j=1}^{\infty} =

%       {\theta_j(x)}_{j=1}^{\infty}
%       U
%       {\psi_j(x) * u}_{j=1}^{\infty}
%       U
%       {u^T R u}
%
%   TRUNCATION:
%   {\theta_j(x)}_{j=1}^{N_3}
%   {\psi_j(x)}_{j=1}^{N_2}
%   N_H = N_2 + N_3
%
%   \Phi(x) = (\phi_1(x), ... , \phi_{N1}(x))
%   \Sigma(x,u) = [     \psi_1(x) u
%                       ...
%                       \psi_{N2}(x) u
%                       \theta_1(x)
%                       ...
%                       \theta_{N3}(x))
%                       u^T R u             ]
%
%   WEIGHTS:
%   w_s \in R^{N_2}
%   v_s \in R^{N_3}
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

% Keeps track of if the algorithm is in the learning phase or the final
% phase
global is_learning;


% Critic NN parameters
global c_s;

% Portion of Hamiltonian NN weight vector corresponding to policy
% approximation
global w_s;

% Portion of Hamiltonian NN weight vector NOT corresponding to policy
% approximation
global v_s;

% Total Hamiltonian weight vector
global wv_s;

% Initial policy u_0(x) (optional)
global u_0;

% Needed for VI weight tuning
global K_phi
global K_sigma;
global ISigmadPhidx;
global ISigmaQ;
global ISigmaR;
global tvec;
global xmat;
global tf;
global x0;


% ***********************
%       
% OPERATION MODE CONTROLS
%

global r1_R0;
global include_in_Sigma;
global int_H_ode45_1_man_0;


% *************************************************************************
% 
% UNPACK ALGORITHM SETTINGS/PARAMETERS
% 
% *************************************************************************


preset = alg_settings.preset;               
sys = alg_settings.sys;         % System

% Get system dimensions
n = sys.n;          % Order of system
m = sys.m;          % Number of inputs

% State, control penalties
Q = alg_settings.Q;
R = alg_settings.R;

% Probing noise
noise = alg_settings.noise;

% Weight differential equation learning time s_f
sf = alg_settings.sf;

% Length of learning window [0, t_f] (sec).
tf = alg_settings.tf;

% Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
% (sec).
tsim = alg_settings.tsim; 

% Max stepsize for ode45
maxstep = alg_settings.maxstep;


% Last function in Hamiltonian basis is r(x,u) (=1) or u^T R u (=0)
r1_R0 = alg_settings.r1_R0;
        
% Whether to include the last Hamiltonian activation function explicitly in
% the basis
include_in_Sigma = alg_settings.include_in_Sigma;

% Whether to integrate Hamiltonian function (13) via ode45 (=1) or manually
% (=0)
int_H_ode45_1_man_0 = alg_settings.int_H_ode45_1_man_0;

% ***********************
%       
% BASIS SETTINGS
%
% See description for details.
%

% Basis struct
basis = alg_settings.basis;

% Dimension of critic NN
N1 = basis.Phi.N;

% Dimensions of Hamiltonian NN
N2 = basis.Psi.N;            % Functions associated with actor network
N3 = basis.Theta.N;          % Remaining functions not associated w/ actor
N_H = N2 + N3;               % Total dimension


% Add on running cost explicitly to basis
if include_in_Sigma
    N_H = N_H + 1;
end


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
% following types of dynamic variables will be necessary for simulation:
%
%   System states x
%   Integrals associated with basis functions
%

% Initial conditions
x0 = alg_settings.x0;

% Initial condition for simulation
% See odefunct() for a description of the state partition
x0_sim = [  x0
            zeros(N1^2, 1)
            zeros(N_H^2, 1)
            zeros(N_H * N1, 1)
            zeros(N_H, 1)
            zeros(N_H, 1)       ];



% ***********************
%       
% MISCELLANEOUS VARIABLES
%

% Initialize critic NN weight vector
c_s = alg_settings.c_0;     % Holds critic NN weights of current iteration

% Initialize Hamiltonian NN weight vector
wv_s = zeros(N_H,1);     % Holds Hamiltonian NN weights of current iteration   

% Set learning flag
is_learning = 1;


% Initial policy u_0(x) (optional)
u_0 = alg_settings.u_0;

% Initialize vector to hold condition number of K_{\phi}(t_f),
% K_{\sigma}(t_f)
cond_A_vec = zeros(2,1);



% *************************************************************************
% 
% DATA STORAGE
% 
% *************************************************************************


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

% Stores critic NN weights c_s at each iteration of the VI algorithm
% (iteration indexed by row)
c_mat = c_s';

% Stores Hamiltonian NN weights [w_s^T v_s]^T at each iteration of the VI
% algorithm (iteration indexed by row)
wv_mat = wv_s';

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
% LEARNING PHASE: COLLECT ONLINE DATA, EVALUATE BASIS FUNCTION INTEGRALS
%
% Here, the system will be simulated with the initial control/noise input
% over the interval [0, t_f]. In addition, the following integrals will be
% performed:
%
%   IPhiPhi         R^{N_1^2}
%   ISigmaSigma     R^{N_H^2}
%   ISigmadPhidx    R^{N_H * N_1}
%   ISigmaQ         R^{N_H}
%   ISigmaR         R^{N_H}
%
% See the odefunct() for a description of these matrices.
% 
% *************************************************************************
% *************************************************************************


% *************************************************************************
%
% COLLECT SIMULATION DATA
%  
% *************************************************************************


% *************************************************************************
%       
% RUN SIMULATION
%

tspan = [0 tf];

% Set max step size for ode45
ode_opts = odeset('MaxStep',maxstep);

[t, x] = ode45(@odefunct, tspan, x0_sim, ode_opts);

% % DEBUGGING: Plot state trajectory
% figure(100)
% hold on
% plot(t,x(:,1:n));
% title('State Trajectory')
% grid on
% xlabel('Time (sec)')
% ylabel('x(t)')

% DEBUGGING: Print IC x(0)
x0

% *************************************************************************
%       
% EXTRACT DYNAMIC VARIABLE DATA
%         

% % Extract dynamic state variable data at beginning of simulation
% x0 = (x(1,:))';

% Extract dynamic state variable data at end of simulation
x1 = (x(size(x,1),:))';

% Unpack state data at end of simulation. See odefunct() for a description
% of the state partition.

% State dynamic variables
xs1 = x1(1:n);       

% Extract dynamic variables associated with integrals of basis functions
% (vectors)
ind0 = n;
len = N1^2;
IPhiPhi = x1(ind0+1:ind0+len);          % R^{N_1^2}
ind0 = ind0 + len;
len = N_H^2;
ISigmaSigma = x1(ind0+1:ind0+len);      % R^{N_H^2}
ind0 = ind0 + len;
len = N_H * N1;
ISigmadPhidx = x1(ind0+1:ind0+len);     % R^{N_H * N_1}
ind0 = ind0 + len;
len = N_H;
ISigmaQ = x1(ind0+1:ind0+len);          % R^{N_H}
ind0 = ind0 + len;
len = N_H;
ISigmaR = x1(ind0+1:ind0+len);          % R^{N_H}

% IPhiPhi = x1(n+1:n+N1^2);                                % R^{N_1^2}
% ISigmaSigma = x1(n+N1^2+1:n+N1^2+N_H^2);                 % R^{N_H^2}
% ISigmadPhidx = x1(n+N1^2+N_H^2+1:n+N1^2+N_H^2+N1*N_H);   % R^{N_H * N_1}
% ISigmaQ = x1(n+N1^2+N_H^2+N1*N_H+1:n+N1^2+N_H^2+N1*N_H+N_H);   % R^{N_H}

% Reshape dynamic variables associated with integrals of basis functions
% into matrices of appropriate dimension
%
% K_\phi(t_f) \in R^{N_1 x N_1} (cf. Sec. IV. A.)
K_phi = reshape(IPhiPhi, [N1, N1]);
% inv_K_phi = inv(K_phi);                 % K_\phi^{-1}(t_f)

% K_\sigma(t_f) \in R^{N_H x N_H} (cf. Sec. IV. A.)
K_sigma = reshape(ISigmaSigma, [N_H, N_H]);
% inv_K_sigma = inv(K_sigma);                 % K_\sigma^{-1}(t_f)

% \int_{0}^{t_f} \Sigma(x,u) [\nabla\Phi(x) \dot{x}]^T dt \in R^{N_H x N_1}
% Used for updating Hamiltonian NN weight c_s (cf. Sec. IV. D.)
ISigmadPhidx = reshape(ISigmadPhidx, [N_H, N1]);

% % K_\sigma^{-1}(t_f) * \int_{0}^{t_f} \Sigma(x,u) [\nabla\Phi(x) \dot{x}]^T
% % dt is used for the Hamiltonian NN weight update. Define it here so
% % the matrix inversion doesn't have to be calculated repeatedly. NOTE:
% % A^{-1} b = A\b in MATLAB. The latter is quicker and more accurate.
% %
% inv_K_sigma_ISigmadPhidx = K_sigma \ ISigmadPhidx;          % R^{N_H x N_1}

% % % DEBUGGING
% % inv_K_psi_IPsidPhidx2 = inv(K_psi) * IPsidPhidx;
% % inv_K_psi = inv(K_psi);

% ***********************
%       
% DATA STORAGE -- SYSTEM STATE DATA
% 

% Store time data
tvec = t;
len_tvec_learn = size(tvec,1);   % Length of time vector for learning phase

% Store system state data
xmat = x(:,1:n);


% Store condition number data
cond_A_vec(1) = cond(K_phi);
cond_A_vec(2) = cond(K_sigma);

% DEBUGGING: Check problem conditioning
disp(['Condition Number of K_{\phi}(t_f) for Least Squares:         '...
            num2str(cond(K_phi), 4)])
disp(['Condition Number of K_{\sigma}(t_f) for Least Squares:       '...
            num2str(cond(K_sigma), 4)])


%%
% *************************************************************************
% *************************************************************************
%
% VI PHASE: CALCULATE CRITIC NN WEIGHTS c_s, HAMILTONIAN NN WEIGHTS w_s,
% v_s               
% 
% *************************************************************************
% *************************************************************************


% *************************************************************************
%       
% RUN SIMULATION
%

tspan = [0 sf];

% % Set max step size for ode45
% ode_opts = odeset('MaxStep',maxstep);

[svec, c_mat] = ode45(@odefunct_vi, tspan, c_s);

% *************************************************************************
%       
% CALCULATE HAMILTONIAN WEIGHTS w_s, v_s BASED ON CRITIC WEIGHTS c_s
%

% Length of VI simulation
lensvec = size(svec, 1);

% Data storage
w_mat = zeros(lensvec, N2);
v_mat = zeros(lensvec, N3);

% Calculate the weights
for i = 1:lensvec
    
    % Critic weight c(s)
    cs = c_mat(i,:)';
    
    % Hamiltonian weights update  
    wvs = update_H_weights(cs);    
    
    % Extract both components of the weight vector
    %
    %   wv_s = [    w_s
    %               v_s   ]
    
    w_mat(i,:) = wvs(1:N2);               % R^{N_2}
    v_mat(i,:) = wvs(N2+1:N2+N3);           % R^{N_3}

end



% DEBUGGING: Final critic NN params
c_s

% DEBUGGING: Final Hamiltonian NN params w_s associated with the actor
w_s

% DEBUGGING: Final Hamiltonian NN params v_s NOT associated with the actor
v_s

% DEBUGGING: For case that u^T R u is included in the basis, print its
% associated final weight
if include_in_Sigma
   runningcost_weight = wv_s(N_H) 
end


%%
% *************************************************************************
% *************************************************************************
%
% FINAL PHASE: APPLY APPROXIMATE OPTIMAL CONTROL
%
%   \mu(x, wv(s_f)) ~ \mu(x, wv_{i*})
%                               
%
% (cf. Alg. 1) where i* denotes the terminating iteration of the VI
% algorithm.
% 
% *************************************************************************
% *************************************************************************

if norm(c_s) < 1e3  % Only simulate if the weights haven't diverged

    % Clear learning flag
    is_learning = 0;

    % ***********************
    %       
    % RUN SIMULATION
    %

    % Initial condition
    x0_sim = xs1;

    % Time span for simulation
    tspan = [tf, tf + tsim];

    % Run simulation
    [t, x] = ode45(@odefunct_final, tspan, x0_sim);


    % ***********************
    %       
    % STORE DATA
    %

    % Store time data
    tvec = [    tvec
                t       ];

    % Store system state data
    xmat = [    xmat
                x       ];
            
else
    
    % Else, the weights have diverged. Warn user
    
    beep
    
    disp('***************************************************************')
    disp('***************************************************************')
    disp('*')
    disp('*')
    disp(' *** WARNING: VI ALGORITHM -- WEIGHTS DIVERGED ***')
    disp('*')
    disp('*')
    disp('***************************************************************')
    disp('***************************************************************')

end
            
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

% Simulation data
out_data.tvec = tvec;
out_data.xmat = xmat;

out_data.lentvec_learn = len_tvec_learn;

% Weights data
out_data.svec = svec;
out_data.c_mat = c_mat;
out_data.w_mat = w_mat;
out_data.v_mat = v_mat;

% Condition number data
out_data.cond_A_vec = cond_A_vec;


% *************************************************************************
%
% CONTROL SIGNAL
% 
% *************************************************************************

% Initialize empty matrix
umat = zeros(size(tvec,1), m);

% Calculate control
for k = 1:size(tvec,1)
    
    % Extract time
    t = tvec(k);
    
    % Extract state
    xs = xmat(k,:)';
    
    % Check if in learning phase and set flag appropriately
    if k <= len_tvec_learn
        is_learning = 1;
    else
        is_learning = 0;
    end    
    
    % Evaluate control 
    u = uxt_alg(xs, t);
                
    % Store control
    umat(k,:) = u';

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
% CALCULATE DYNAMICS FOR LEARNING PORTION OF ALGORITHM
%
% NOTE: State partitioning for simulation
%
% In order to perform the VI algorithm, the following two types of dynamic
% variables will be necessary for simulation:
%
%   System states
%   Integrals associated with activation functions
% 
% *************************************************************************
%
% DYNAMIC VARIABLES ASSOCIATED WITH INTEGRALS OF ACTIVATION FUNCTIONS
%
% Here, the system will be simulated with the initial control/noise input
% over the interval [0, t_f]. In addition to the system state variables,
% the following integrals will be performed:
%
%   IPhiPhi         R^{N_1^2}
%   ISigmaSigma     R^{N_H^2}
%   ISigmadPhidx    R^{N_H * N_1}
%   ISigmaQ         R^{N_H}
%   ISigmaR         R^{N_H}
%
% A description of each of these matrices can be found below:
%
%                   IPhiPhi     R^{N_1^2}
%
% This integral is given by,
%
%   IPhiPhi = \int_{0}^{t_f} kron(\Phi(x), \Phi(x)) dt
%
% where kron(a,b) denotes the Kronecker tensor product of two vectors a, b.
% IPhiPhi is used to calculate the matrix (cf. Sec. IV. A.),
%
%   K_\phi(t_f) \in R^{N_1 x N_1},
%   K_\phi(t_f) = \int_{0}^{t_f} \Phi(x) \Phi^T(x) dt 
%
%                   ISigmaSigma     R^{N_H^2}
%
% This integral is given by,
%
%   ISigmaSigma = \int_{0}^{t_f} kron(\Sigma(x,u), \Sigma(x,u)) dt
%
% where kron(a,b) denotes the Kronecker tensor product of two vectors a, b.
% ISigmaSigma is used to calculate the matrix (cf. Sec. IV. A.),
%
%   K_\sigma(t_f) \in R^{N_H x N_H},
%   K_\sigma(t_f) = \int_{0}^{t_f} \Sigma(x,u) \Sigma^T(x,u) dt 
%   
%                   ISigmadPhidx  R^{N_H * N_1}
%
% This integral is given by,
%
%   IPsidPhidx = vec( \Sigma(x,u) [\nabla\Phi(x) \dot{x}]^T dt )
%
% where vec(A) denotes the vectorization of the matrix A. This has to be
% performed because ode45 can only integrate a vector of states, not a
% matrix of states. ISigmadPhidx is used to calculate the matrix (cf. Sec.
% IV. A.),
%
%   \int_{0}^{t_f} \Sigma(x,u) [\nabla\Phi(x) \dot{x}]^T dt 
%                   \in R^{N_H x N_1}
%
% which, it can be checked, given c \in R^{N_1}, satisfies the identity:
%
%   [ \int_{0}^{t_f} \Sigma(x,u) [\nabla\Phi(x) \dot{x}]^T dt ] * c
%   =
%   \int_{0}^{t_f} \Sigma(x,u) {d \hat{V}_{N1}(x, w)}.
%
% Where \nabla\Phi(x) \in R^{N_1 x n} is the Jacobian matrix of the critic
% activation functions \Phi(x). Given that the control penalty U^T R u is a
% member of the Hamiltonian NN basis, the right-hand-side integral is
% exactly the integral involved in the Hamiltonian NN weight update
% (cf. Sec. IV. A.).
%
%                   ISigmaQ     R^{N_H}
%
% This integral is given by,
%
%   ISigmaQ = \int_{0}^{t_f} \Sigma(x,u) * Q(x) dt
%
%                   ISigmaR     R^{N_H}
%
% This integral is given by,
%
%   ISigmaR = \int_{0}^{t_f} \Sigma(x,u) * u^T R u dt
%
% *************************************************************************
%
% FINAL STATE PARTITION:
%
% x_sim =   [   x                   R^n
%               IPhiPhi             R^{N_1^2}
%               ISigmaSigma         R^{N_H^2}
%               ISigmadPhidx        R^{N_H * N_1}
%               ISigmaQ             R^{N_H} 
%               ISigmaR             R^{N_H}         ]
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

function xdot = odefunct(t, x)

% Global variables
global sys;
global Q;
global R;
global basis;

global r1_R0;
global include_in_Sigma;

% Get system dimensions
n = sys.n;          % Order of system
% m = sys.m;          % Number of inputs
    
 
% ***********************
%       
% SYSTEM DYNAMICS
%

% Extract system variables
xs = x(1:n);

% Evaluate system drift dynamics
fx = eval_f(xs, sys);

% Evaluate input gain matrix
gx = eval_g(xs, sys);

% Evaluate control signal u(t)
u = uxt_alg(xs, t);

% Evaluate state derivative \dot{x}
dxs = fx + gx * u;

% % % ***********************
% % %       
% % % RUNNING COST r(x,u)
% % %
% % 
% % % Evaluate state penalty Q(x)
% % Qx = eval_Q(xs, Q);
% % 
% % % Evaluate control penalty
% % Rxu = eval_R(xs, u, R);
% % 
% % % Evaluate running cost r(x,u)
% % rxu = Qx + Rxu;


% ***********************
%       
% BASIS FUNCTIONS
%

% Dimension of critic NN
N1 = basis.Phi.N;

% Dimensions of Hamiltonian NN
N2 = basis.Psi.N;            % Functions associated with actor network
N3 = basis.Theta.N;          % Remaining functions not associated w/ actor
N_H = N2 + N3;               % Total dimension

% Add on running cost explicitly to basis
if include_in_Sigma
    N_H = N_H + 1;
end

% Evaluate basis functions
[Phix, dPhix] = eval_phi(xs, basis.Phi);

% \Sigma(x,u)
% BEFORE multiplying first N_2 functions on the right by u 
Psix = eval_phi(xs, basis.Psi);         % First N_2 components
Thetax = eval_phi(xs, basis.Theta);     % Last N_3 components

% Multiply factor of u in first N_2 elements
Psixu = Psix * u;

% Combine activation functions into one vector
Sigmaxu = [     Psixu 
                Thetax  ];

% Evaluate state penalty Q(x)
Qx = eval_Q(xs, Q);            
            
% Evaluate control penalty
Rxu = eval_R(xs, u, R);  

% Evaluate running cost r(x,u)
rxu = Qx + Rxu;  

            
% Add on running cost explicitly to basis
if include_in_Sigma
    if r1_R0
        Sigmaxu = [Sigmaxu ; rxu]; 
    else
        Sigmaxu = [Sigmaxu ; Rxu]; 
    end   
end


% ***********************
%       
% DYNAMIC VARIABLES ASSOCIATED WITH BASIS FUNCTION INTEGRALS
%

% For K_phi = \int_{0}^{tf} Phi(x) Phi^T(x) dt
% cf. Sec. IV. A.
dIPhiPhi = kron(Phix, Phix);

% For K_psi = \int_{0}^{tf} \Sigma(x,u) \Sigma^T(x,u) dt
% cf. Sec. IV. A.
dISigmaSigma = kron(Sigmaxu, Sigmaxu);

% For \int_{0}^{tf} \Sigma(x,u) [\nabla \Phi(x) \dot{x}]^T dt
% cf. Sec. IV. D.
% Original shape (N_H x N1)
dISigmadPhidx = Sigmaxu * (dPhix * dxs)';
% Reshaped into an N_H*N1-dimensional vector
dISigmadPhidx = reshape(dISigmadPhidx, N_H*N1, 1);

% For ISigmaQ = \int_{0}^{t_f} \Sigma(x,u) * Q(x) dt
dISigmaQ = Sigmaxu * Qx;

% For ISigmaR = \int_{0}^{t_f} \Sigma(x,u) * u^T R u dt
dISigmaR = Sigmaxu * Rxu;


% ***********************
%       
% PACKAGE STATE DERIVATIVE OUTPUT
%
% x_sim =   [   x                   R^n
%               IPhiPhi             R^{N_1^2}
%               ISigmaSigma         R^{N_H^2}
%               ISigmadPhidx        R^{N_H * N_1}
%               ISigmaQ             R^{N_H}         
%               ISigmaR             R^{N_H}         ]
%
xdot = [    dxs
            dIPhiPhi
            dISigmaSigma
            dISigmadPhidx
            dISigmaQ
            dISigmaR        ];
        
        
%%
% *************************************************************************
% *************************************************************************
%
% CALCULATE WEIGHT DYNAMICS FOR VI PHASE: CRITIC NN WEIGHTS c_s,
% HAMILTONIAN NN WEIGHTS w_s, v_s
%
% We shall use ode45 for weight simulation, NOT the forward-Euler
% approximation of the updates (13) and (14) found at the beginning of Sec.
% IV. D.
%
% The Hamiltonian NN weight update is as follows:
%
%   wv_s = K_\sigma^{-1}(t_f) \int_{0}^{t_f} \Sigma(x,u) * 
%                           (d \hat{V}_{N1}(x, c_s) + r) dt
%
% where wv_s = [w_s^T v_s^T]^T is shorthand
%
% The critic NN weight update is as follows:
%
% d/ds c_s = K_\phi^{-1}(t_f) \int_{0}^{t_f} \Phi(x) *
%         \hat{H}_{N_H}(x, \mu_{N_H}(x, wv_s), wv_s) dt
%
% The integral in the above will have to be done manually with the
% pre-existing generated state dynamic data, given that the Hamiltonian
% \hat{H}_{N_H}(x, \mu_{N_H}(x, wv_s), wv_s) is nonlinear in the weights
% wv_s = [v_s^T w_s^T]^T. This is the reason why the associated integral
% could not be performed in the learning phase, as was done with the
% Hamiltonian weights update integral. See odefunct() for a more detailed
% explanation.
%
% *************************************************************************
%
% FINAL STATE PARTITION:
%
% x_sim =   [   c_s                 R^{N_1}  ]
%                               
% 
% *************************************************************************
% *************************************************************************        

function xdot = odefunct_vi(t, x)

% Global variables
global K_phi
% global K_sigma;
% global ISigmadPhidx;
% global ISigmaQ;
% global ISigmaR;

global r1_R0;
global include_in_Sigma;
global int_H_ode45_1_man_0;

global tvec;
global xmat;

global R;
global Q;
global basis;
global tf;
global x0;

global c_s;
global w_s;
global v_s;
global wv_s;

%     % DEBUGGING: Print iteration count
%     disp(['s = ' num2str(t)])

% *********************************************************************
%       
% EXTRACT STATE, BASIS PARAMETERS
%

% Extract the critic parameters
c_s = x;


% ***********************
%       
% BASIS FUNCTIONS
%

% Dimension of critic NN
N1 = basis.Phi.N;

% Dimensions of Hamiltonian NN
N2 = basis.Psi.N;            % Functions associated with actor network
N3 = basis.Theta.N;         % Remaining functions not associated w/ actor
N_H = N2 + N3;               % Total dimension

% Add on running cost explicitly to basis
if include_in_Sigma
    N_H = N_H + 1;
end


% *********************************************************************
%       
% HAMILTONIAN NN WEIGHT w_s, v_s UPDATE
%
% cf. Sec. IV. D.
%

% Hamiltonian weights update  
wv_s = update_H_weights(c_s);

% Extract both components of the weight vector. Recall the weight
% vector is partitioned according to the number of basis functions N_2
% chosen in the actor basis {\psi_i(x) * u}_{j=1}^{N_2} and the number
% of basis functions N_3 chosen in the other basis not pertaining to
% the actor {\theta_i(x)}_{j=1}^{N_3}. In sum
%
%   wv_s = [    w_s
%               v_s   ]
%

w_s = wv_s(1:N2);               % R^{N_2}
v_s = wv_s(N2+1:N2+N3);           % R^{N_3}



% *********************************************************************
%       
% CRITIC NN WEIGHT c_s UPDATE
%
% cf. Sec. IV. D.
%


if int_H_ode45_1_man_0
    
    % Perform integral 
    % \int_{0}^{t_f} \Phi(x) * \hat{H}_{N_H}(x, \mu_{N2}(x, wv_s), wv_s) dt
    % via ode45

    tspan = [0 tf];
    x0_sim = [  x0
                zeros(N1, 1)    ];      
    [~, x_IPhiH] = ode45(@odefunct_H, tspan, x0_sim);

    % % DEBUGGING: Plot state trajectory
    % [t, x_IPhiH] = ode45(@odefunct_H, tspan, x0_sim);
    % figure(101)
    % hold on
    % plot(t,x_IPhiH(:,1:end-N1));
    % title('State Trajectory -- Hamiltonian Integral')
    % grid on
    % xlabel('Time (sec)')
    % ylabel('x(t)')


    % Extract integral
    IPhiH = x_IPhiH(end, end-N1+1:end)';
    

else
    
    % Perform integral 
    % \int_{0}^{t_f} \Phi(x) * \hat{H}_{N_H}(x, \mu_{N2}(x, wv_s), wv_s) dt
    % Manually with pre-existing state trajectory data

    IPhiH = zeros(N1,1);  % Integration variable
    
    % ***********************
    %       
    % INTEGRATION LOOP
    %   
    for k = 1:size(tvec,1)-1

        % Calculate time differential
        dt = tvec(k+1) - tvec(k);

        % Get state vector at current time value
        xs = xmat(k,:)';

        % ***********************
        %       
        % BASIS FUNCTIONS
        %        

        % \Phi(x)
        Phix = eval_phi(xs, basis.Phi);


        % \Sigma(x,u)
        % BEFORE multiplying first N_2 functions on the right by u 
        Psix = eval_phi(xs, basis.Psi);         % First N_2 components
        Thetax = eval_phi(xs, basis.Theta);     % Last N_3 components

        % ***********************
        %       
        % EVALUATE POLICY \mu(x, wv_k)
        %    
        mux = - (1 / 2) * inv(R) * (w_s' * Psix)';


        % Multiply factor of \mu(x, wv_k) in first N_2 elements
        Psixmu = Psix * mux;

        % Combine activation functions into one vector
        Sigmaxmu = [    Psixmu 
                        Thetax  ];                  

        % ***********************
        %       
        % RUNNING COST r(x,\mu_{N_H}(x, wv_k))
        %

        % Evaluate state penalty Q(x)
        Qx = eval_Q(xs, Q);

        % Evaluate control penalty
        Rxmu = eval_R(xs, mux, R);

        % Evaluate running cost r(x,\mu_{N_H}(x, wv_k))
        rxmu = Qx + Rxmu;     

        % Add on running cost explicitly to basis
        if include_in_Sigma
            if r1_R0
                Sigmaxmu = [Sigmaxmu ; rxmu];
            else
                Sigmaxmu = [Sigmaxmu ; Rxmu];
            end   
        end

        % ***********************
        %       
        % EVALUATE HAMILTONIAN APPROXIMATION 
        %
        % \hat{H}(x, \mu(x, wv_k), wv_k)
        %
        %       = wv_k^T * \Sigma(x,\mu(x, wv_k)) + R(\mu(x, wv_k))
        %

        if include_in_Sigma
            H_x_mu = wv_s' * Sigmaxmu;
        else
            if r1_R0
                H_x_mu = wv_s' * Sigmaxmu + rxmu;
            else
                H_x_mu = wv_s' * Sigmaxmu + Rxmu;
            end
        end


        % ***********************
        %       
        % INCREMENT INTEGRAND
        %    

        IPhiH = IPhiH + Phix * H_x_mu * dt;

    end

end



% ***********************
%       
% WEIGHT DERIVATIVE
%      

dc = K_phi \ IPhiH;

 
% ***********************
%       
% PACKAGE STATE DERIVATIVE OUTPUT
%
% x_sim =   [   c_s                 R^{N_2}  ]
%                
   
xdot = [    dc  ];



%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INTEGRATE HAMILTONIAN FOR VI PORTION OF ALGORITHM
%
% This is for performing the integral,
%
% \int_{0}^{t_f} \Phi(x) * \hat{H}_{N_H}(x, \mu_{N2}(x, wv_s), wv_s) dt
%
% which is associated with the critic differential weight update
%
%   d / ds {c(s)}.
%
% The state variables x as well as the Hamiltonain integral above will need
% to be simulated as dynamic variables.
%
% NOTE: The state trajectry generated {x(t)}_{t=0}^{t_f} is the same as
% that generated by the initial learning phase. This separate integral was
% found necessary for numerics purposes.
%
% *************************************************************************
%
% FINAL STATE PARTITION:
%
% x_sim =   [   x           R^n     
%               IPhiH       R^{N_1}     ]
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

function xdot = odefunct_H(t, x)


% Global variables
global sys;
global R;
global Q;
global basis;

global r1_R0;
global include_in_Sigma;


% global c_s;
global w_s;
% global v_s;
global wv_s;


% Get system dimensions
n = sys.n;          % Order of system
% m = sys.m;          % Number of inputs

% ***********************
%       
% SYSTEM DYNAMICS
%

% Extract system variables
xs = x(1:n);

% Evaluate system drift dynamics
fx = eval_f(xs, sys);

% Evaluate input gain matrix
gx = eval_g(xs, sys); 

% Evaluate control u(t)
u = uxt_alg(xs, t);

% Evaluate state derivative \dot{x}
dxs = fx + gx * u;


% ***********************
%       
% HAMILTONIAN INTEGRAL
%        

% \Phi(x)
Phix = eval_phi(xs, basis.Phi);


% \Sigma(x,u)
% BEFORE multiplying first N_2 functions on the right by u 
Psix = eval_phi(xs, basis.Psi);         % First N_2 components
Thetax = eval_phi(xs, basis.Theta);     % Last N_3 components

% ***********************
%       
% EVALUATE POLICY \mu(x, wv_k)
%    
mux = - (1 / 2) * inv(R) * (w_s' * Psix)';


% Multiply factor of \mu(x, wv_k) in first N_2 elements
Psixmu = Psix * mux;

% Combine activation functions into one vector
Sigmaxmu = [    Psixmu 
                Thetax  ];                  

% ***********************
%       
% RUNNING COST r(x,\mu_{N_H}(x, wv_k))
%

% Evaluate state penalty Q(x)
Qx = eval_Q(xs, Q);

% Evaluate control penalty
Rxmu = eval_R(xs, mux, R);

% Evaluate running cost r(x,\mu_{N_H}(x, wv_k))
rxmu = Qx + Rxmu;     

% Add on running cost explicitly to basis
if include_in_Sigma
    if r1_R0
        Sigmaxmu = [Sigmaxmu ; rxmu];
    else
        Sigmaxmu = [Sigmaxmu ; Rxmu];
    end   
end

% ***********************
%       
% EVALUATE HAMILTONIAN APPROXIMATION 
%
% \hat{H}(x, \mu(x, wv_k), wv_k)
%
%       = wv_k^T * \Sigma(x,\mu(x, wv_k)) + R(\mu(x, wv_k))
%

if include_in_Sigma
    H_x_mu = wv_s' * Sigmaxmu;
else
    if r1_R0
        H_x_mu = wv_s' * Sigmaxmu + rxmu;
    else
        H_x_mu = wv_s' * Sigmaxmu + Rxmu;
    end
end


% ***********************
%       
% HAMILTONIAN DIFFERENTIAL
%    

dIPhiH = Phix * H_x_mu;


% ***********************
%       
% PACKAGE STATE DERIVATIVE OUTPUT
%
% x_sim =   [   x           R^n     
%               IPhiH       R^{N_1}     ]
%
xdot = [    dxs
            dIPhiH  ];





%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CALCULATE DYNAMICS FOR FINAL PORTION OF ALGORITHM
%
% Note now that the dynamic variables associated with the activation
% function integrals no longer need to be simulated, as learning has
% concluded. The state is then simply the dynamic system state.
%
% *************************************************************************
%
% FINAL STATE PARTITION:
%
% x_sim =   [   x           R^n     ]
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

function xdot = odefunct_final(t, x)

% Global variables
global sys;

% Get system dimensions
n = sys.n;          % Order of system
% m = sys.m;          % Number of inputs

% ***********************
%       
% SYSTEM DYNAMICS
%

% Extract system variables
xs = x(1:n);

% Evaluate system drift dynamics
fx = eval_f(xs, sys);

% Evaluate input gain matrix
gx = eval_g(xs, sys); 

% Evaluate policy \mu_{N_2}(x, c_s)
u = uxt_alg(xs, t);

% Evaluate state derivative \dot{x}
dxs = fx + gx * u;


% ***********************
%       
% PACKAGE STATE DERIVATIVE OUTPUT
%
% x_sim =   [   x           R^n     ]
%
xdot = [    dxs    ];


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


function u = uxt_alg(x, t)

% Global variables
% global preset;
global sys;
global R;
global w_s;
global basis;
global noise;
global u_0;

global is_learning;

% Get system dimensions
% n = sys.n;          % Order of system
m = sys.m;


if is_learning
    
    % *********************************************************************
    %
    % LEARNING PHASE -- APPLY PROBING NOISE
    %
	% *********************************************************************

    
    % Evaluate noise e(t)
    et = eval_noise(t, m, noise);

    % Evaluate control signal u(t)
    u = eval_u(x, u_0) + et;    
    
    
else

    % *********************************************************************
    %
    % POST-LEARNING
    %
	% *********************************************************************    

    % ***********************
    %       
    % BASIS FUNCTIONS
    %        

    % \Sigma(x,u)
    % BEFORE multiplying first N_2 functions on the right by u 
    Psix = eval_phi(x, basis.Psi);         % First N_2 components
%     Thetax = eval_phi(xs, basis.Theta);     % Last N_3 components

%     % Combine activation functions into one vector
%     Sigmaxu = [     Psixu 
%                     Thetax  ];

    % ***********************
    %       
    % EVALUATE POLICY \mu_{N_H}(x, wv_s)
    %   
    u = - (1 / 2) * inv(R) * (w_s' * Psix)';


end


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% HAMILTONIAN WEIGHT UPDATE
%
% See (14), T. Bian and Z.-P. Jiang.
%
% The Hamiltonian NN weight update is as follows:
%
%   wv_s = K_\sigma^{-1}(t_f) \int_{0}^{t_f} \Sigma(x,u) * 
%                           (d \hat{V}_{N1}(x, c_s) + r) dt
%
% where wv_s = [w_s^T v_s^T]^T is shorthand.
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************


function wv_s = update_H_weights(cs)


% Global variables
% global K_phi
global K_sigma;
global ISigmadPhidx;
global ISigmaQ;
global ISigmaR;

global r1_R0;
global include_in_Sigma;

if include_in_Sigma
    
    wv_s = K_sigma \ (ISigmadPhidx * cs + ISigmaQ + ISigmaR);
    
else
    
    if r1_R0
        
        wv_s = K_sigma \ (ISigmadPhidx * cs);
        
    else
        
        wv_s = K_sigma \ (ISigmadPhidx * cs + ISigmaQ);
%         wv_s = K_sigma \ (ISigmadPhidx * cs);
%         wv_s = K_sigma \ (ISigmadPhidx * cs + ISigmaQ + ISigmaR);
        
    end
    
end



