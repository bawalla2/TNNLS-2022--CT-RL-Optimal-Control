function out_data = alg_irl(alg_settings)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INTEGRAL REINFORCEMENT LEARNING (IRL) ALGORITHM
%
% Brent Wallace  
%
% 2021-11-06
%
% This program implements the IRL algorithm presented in,
%
%   D. Vrabie and F.L. Lewis. "Neural network approach to continuous-time
%   direct adaptive optimal control for partially unknown nonlinear
%   systems." Neural Networks, 22:237-246, 2009.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% out_data = alg_irl(alg_settings)
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
%   basis                   (Struct) contains activation function basis
%                           parameters. Has the following fields:
%       .tag                (String) tag of the desired activation function
%                           basis to use (see eval_phi.m for options).
%       .N                  (Integer) the integer "N," or the basis
%                           dimension.
%   noise                   (Struct) contains info for probing noise. NOTE:
%                           Not required for this algorithm. Has the
%                           following fields:
%       .tag                (String) tag of specific probing noise signal
%                           to be injected (see eval_noise.m for options).
%                           If no noise is desired, simply enter this as
%                           '0'.
%   T                       (Double) integral reinforcement interval length
%                           (sec).
%   istar                   (Integer) number of policy iterations to
%                           execute before terminating the algorithm.
%   num_sims_per_iter       (Integer) Number of simulations to execute per
%                           iteration of the PI algorithm. NOTE: a
%                           simulation consists of 'l'
%                           data samples, each taken over 'T' seconds. So
%                           each iteration has 'num_sims_per_iter' *
%                           'l' total T-second samples
%                           collected.
%   l     (Integer) number of samples to collect per
%                           simulation.     
%   x0mat                   (Matrix) This matrix can be empty or have up to
%                           'istar' * 'num_sims_per_iter' rows, with
%                           n columns. Each time a new simulation begins,
%                           the initial conditions need to be set. These
%                           ICs can be set manually for as many of the
%                           total 'istar' * 'num_sims_per_iter'
%                           simulations run in the algorithm. If x0mat runs
%                           out of ICs, then the ICs will be generated
%                           either randomly at the beginning of each new
%                           simulation, or will be carried over from the
%                           previous simulation's final values (see
%                           variable x0_behavior for details).
%   x0_behavior             (String) Determines the manner in which ICs are
%                           generated at the beginning of each new
%                           simulation, after the algorithm has run out of
%                           ICs specified in x0mat. Has the options:
%       'rand'              Randomly generate new ICs.
%       'cont'              Continue ICs of new simulation as final values
%                           of previous simulation.
%   c_0                      (N-dimensional vector) ICs for critic NN, where
%                           N is the critic NN basis dimension.
%   tsim                    (Double) time window for post-learning
%                           simulation (sec). I.e., if learning happens
%                           over [0, t_f], post-learning happens over [t_f,
%                           t_f + tsim].
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% out_data                  (Struct) algorithm output data. Has the
%                           following fields:
%       .tvec               ('simlength'-dimensional vector) vector of time
%                           indices corresponding to the simulation time
%                           instants over the course of the algorithm
%                           execution.
%       .xmat               ('simlength' x n matrix) Matrix whose row
%                           indexes the time instants specified in .tvec,
%                           and whose n-columns are the state vector at the
%                           respective time instant.
%       .umat               ('simlength' x m matrix) Matrix whose row
%                           indexes the time instants specified in .tvec,
%                           and whose m-columns are the control signal u(t)
%                           at the respective time instant.
%       .tvec_pi            ('istar' - dimensional vector)
%                           Vector of time instants corresponding to the
%                           sample instants of each new iteration of the PI
%                           algorithm.
%       .c_mat               (N x 'istar' matrix) critic NN weights
%                           at each of the time instants of .tvec_pi.
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

global Q;
global R;
global sys;

global basis;
global noise;

global c_i;

% Keeps track of if the algorithm is in the learning phase or the final
% phase
global is_learning;



% *************************************************************************
% 
% UNPACK ALGORITHM SETTINGS/PARAMETERS
% 
% *************************************************************************

T = alg_settings.T;
istar = alg_settings.istar;
num_sims_per_iter = alg_settings.num_sims_per_iter;
l = alg_settings.l;
        
sys = alg_settings.sys;             % System array
n = alg_settings.sys.n;             % System order

Q = alg_settings.Q;
R = alg_settings.R;

basis = alg_settings.basis;         % Basis struct
N = alg_settings.basis.N;           % Basis dimension

% Probing noise
noise = alg_settings.noise; 

% Post-learning simulation length
tsim = alg_settings.tsim;


% *************************************************************************
% 
% ALGORITHM INITIALIZATION
% 
% *************************************************************************

% ***********************
%       
% MISCELLANEOUS VARIABLES
%

% Matrix of ICs
x0mat = alg_settings.x0mat;

% IC generation behavior. See above for details.
x0_behavior = alg_settings.x0_behavior;

% Initial critic NN weights
c_i = alg_settings.c_0;

% Set learning flag
is_learning = 1;

% Initialize simulation counter
simcount = 1;

% Initialize vector to hold condition number at each iteration
cond_A_vec = zeros(istar, 1);


% *************************************************************************
% 
% DATA STORAGE
% 
% *************************************************************************


% Critic NN weights
c_mat = zeros(istar + 1,N);
c_mat(1,:) = c_i';

% Time vector, state trajectory, control signal
tvec = [];
xmat = [];
umat = [];



% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% BEGIN MAIN
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************


for i = 1:istar
% for i = 1:istar - 1          % DEBUGGING

    
    % ***********************
    %       
    % INITIALIZE SAMPLE COUNTER
    %
    % Each iteration of the PI algorithm, 'num_sims_per_iter' *
    % 'l' samples of the state and integral reinforcement
    % are collected for the least squares minimization at the end of the
    % iteration, each data point corresponding to a T-second simulation.
    % This counter is initialized to 1 here and is incremented after each
    % T-second simulation until it will eventually reach the value
    % 'num_sims_per_iter' * 'l'
    %
    samplecount = 1;
    
    % Storage for time data and simulation data at the current PI index i
    tvec_i = [];
    xmat_i = [];
    
    % *********************************************************************
    %
    % OUTER PI LOOP
    %
    % This outer loop corresponds to one iteration of the PI algorithm.
    % 'num_sims_per_iter' simulations are run, and then the critic
    % parameters are updated via least squares.
    %        
    for m = 1:num_sims_per_iter
 
        % ***********************
        %       
        % INITIAL CONDITION FOR SIMULATION
        %
        % Dynamic variables include state x (first n entries) and integral
        % reinforcement ((n+1)-th entry)
        %
        
        % Determine if state IC is to be set manually or by the algorithm
        if simcount <= size(x0mat, 1)
            % Set IC manually
            xs0 = x0mat(simcount,:)';
        else
            % No more manual ICs
            if (samplecount*i == 1) || (strcmp(x0_behavior, 'rand'))
                % Set IC randomly
                xs0 = 2 * (rand(n, 1) - 1/2);    
            else
                % Carry IC over from previous simulation
                xs0 = x1;
            end
        end
        
        % Total IC vector (consists of state + integral reinforcement)
        x0_sim = [  xs0
                    0   ];
        
        
        % *************************************************************
        %
        % INNERMOST LOOP
        %
        % This inner loop is the one in which the individual
        % simulations are run. One simulation will consist of
        % 'l' T-second simulations, each run
        % separately so that data may be collected at the end of the
        % T-second interval.
        %  
        for s = 1:l
            
            % ***********************
            %       
            % RUN SIMULATION
            %
            
            % Time span for current simulation
            tspan = ...
                T * (i - 1) * num_sims_per_iter * l +...
                T * [samplecount - 1, samplecount];

            % Run simulation
%             [t, x] = ode23(@odefunct, [0 T], x0);
            [t, x] = ode45(@odefunct, tspan, x0_sim);
 
%             % DEBUGGING: Plot state trajectory
%             figure(100)
%             plot(t,x(:,1:n));
%             title('State Trajectory')
%             grid on
%             xlabel('Time (sec)')
%             ylabel('x(t)') 
%             hold on
            
%             % DEBUGGING: Plot state trajectory
%             figcount = 200;
%             for xi = 1:n
%                 figure(figcount)
%                 hold on
%                 plot(t,x(:,xi));
%                 title(['State Trajectory x_{' num2str(xi) '}(t)'])
%                 figcount = figcount + 1;
%             end
          
            
            % ***********************
            %       
            % DATA MANAGEMENT
            %         
            
            % Extract state data at end of simulation
            x1 = (x(length(x),(1:end-1)))';
            
            % Extract integral reinforcement at end of simulation
            V(samplecount,1) = x(length(x), end);
            
            % Evaluate the basis function difference
            %
            % phi(x(t + T)) - phi(x(t))
            %
            [phix_tpT, ~] = eval_phi(x1, basis);     % phi(x(t + T))
            [phix_t, ~] = eval_phi(x0_sim(1:n), basis);       % phi(x(t))
            phidiff(samplecount, :) = phix_tpT - phix_t;             
            
            % ***********************
            %       
            % DATA STORAGE
            %              

            tvec_i = [  tvec_i
                        t       ];
                    
            xmat_i = [  xmat_i
                        x(:,1:end-1)    ];
            
            % ***********************
            %       
            % PREPARE NEXT SAMPLE
            % 
            
            % IC for next T-second simulation
            % State is carried over from previous T-second simulation,
            % integral reinforcement is reset
            x0_sim = [ x1 ; 0];            
            
            % Increment sample counter
            samplecount = samplecount + 1;
            
        end
        
        % ***********************
        %       
        % PREPARE NEXT SIMULATION
        %
        
        % Increment simulation counter
        simcount = simcount + 1;
               
    end

    % ***********************
    %       
    % CALCULATE CONTROL SIGNAL APPLIED FOR THIS ITERATION
    %      
    for k = 1:size(tvec_i, 1)
       
        % Get time
        tk = tvec_i(k);
        
        % Get system state
        xs = xmat_i(k,:)';
        
        % Evaluate control signal u(t)
        u = uxt_alg(xs, tk);
        
        % Store control signal
        umat = [    umat
                    u'      ];
        
    end
    
    % ***********************
    %       
    % LEAST SQUARES MINIMIZATION TO DETERMINE CRITIC NN WEIGHTS FOR NEXT
    % ITERATION
    %
    
    c_im1 = c_i;
    
    % Least-squares minimization
    A = - phidiff;
    c_i = A \ V;
    
    % Store new weight
    c_mat(i + 1,:) = c_i';
    
    % Store condition number of matrix involved in least-squares
    % minimization
    cond_A_vec(i) = cond(A);
    
    % DEBUGGING: Show iteration count, norm difference
    disp('*****')
    disp(['i = ' num2str(i) ',     ||c_i - c_{i-1}|| =     ' ...
        num2str(norm(c_i - c_im1))])    
    
    % DEBUGGING: Check problem conditioning
    disp(['Condition Number of "A" for Least Squares:           '...
                num2str(cond(A), 4)])
%     disp(['Condition Number of "A^T A" for Least Squares:       '...
%                 num2str(cond(A'*A), 4)])


    % ***********************
    %       
    % DATA STORAGE
    %              

    tvec = [    tvec
                tvec_i      ];

    xmat = [    xmat
                xmat_i      ];    
    
end

% % DEBUGGING: Plot critic NN params
% figure(101)
% tmp = (0:1:istar) * T * num_sims_per_iter * l;
% plot(tmp, c_mat(:,1:istar + 1)','.-'); 
% title('NN parameters'); 
% xlabel('Time (s)'); 
% lgnd = {'w_1', 'w_2', 'w_3', 'w_4', 'w_5', 'w_6', 'w_7', 'w_8'};
% lgd = legend(lgnd);
% set(lgd, 'Numcolumns', 2);          % Make legend 2 columns
% set(lgd, 'Location', 'Best');       % Put legend in empty spot


% DEBUGGING: ICs
if ~isempty(x0mat)
    x_0 = x0mat(1,:)'
end

% DEBUGGING: Final critic NN params
c_i

%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% POST-LEARNING PHASE: SIMULATE SYSTEM WITH FINAL WEIGHTS
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

% Initial condition. Note now that the integral reinforcement dynamic
% variable no longer needs to be simulated, since the weights are not being
% updated.
x0_sim = x1;

% Total length of learning window [0, t_f]
tf = T * istar * num_sims_per_iter * l;

% Time span for simulation
tspan = [tf, tf + tsim];
% tspan = [tvec(end), tvec(end) + tsim];
        
% Run simulation
[t, x] = ode45(@odefunct_final, tspan, x0_sim);
       
% ***********************
%       
% CALCULATE CONTROL SIGNAL APPLIED POST-LEARNING
%      
for k = 1:size(t, 1)

    % Get time
    tk = t(k);

    % Get system state
    xs = x(k,:)';

    % Evaluate control signal u(t)
    u = uxt_alg(xs, tk);

    % Store control signal
    umat = [    umat
                u'      ];

end


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

% Time, state, control data
out_data.tvec = tvec;
out_data.xmat = xmat;
out_data.umat = umat;

% Time indices corresponding samples of new PI iterations.
tvec_pi = (0:1:istar) * T * num_sims_per_iter * l;

% Weight data
out_data.tvec_pi = tvec_pi;
out_data.c_mat = c_mat;

% Condition number data
out_data.cond_A_vec = cond_A_vec;


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
% CALCULATE DYNAMICS -- LEARNING PHASE
%
% State consists of system state x (n-dimensional), plus the integral
% reinforcement (appended at the (n+1)-th entry).
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

function xdot = odefunct(t, x)

% Global variables
global Q;
global R;
global sys;

% Dynamic state variables
xs = x(1:end-1);

% Evaluate drift dynamics
fx = eval_f(xs, sys);

% Evaluate input gain matrix
gx = eval_g(xs, sys);

% Calculate control signal
u = uxt_alg(xs, t);

% Evaluate state penalty Q(x)
Qx = eval_Q(xs, Q);

% Calculate state derivative
xdot = [    fx + gx * u
            Qx + u' * R * u    ];
        

%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CALCULATE DYNAMICS -- POST-LEARNING PHASE
%
% State consists of system state x (n-dimensional) only.
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

function xdot = odefunct_final(t, x)

% Global variables
global sys;

% Evaluate drift dynamics
fx = eval_f(x, sys);

% Evaluate input gain matrix
gx = eval_g(x, sys);

% Calculate control signal
u = uxt_alg(x, t);

% Calculate state derivative
xdot = fx + gx * u;
        
                
        
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
global sys;
global R;
global basis;
global noise;
global c_i;

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
    
    % Evaluate noise
    et = eval_noise(t, m, noise);
    
    % Evaluate input gain matrix
    gx = eval_g(x, sys);
    
    % Evaluate basis functions and gradient
    [~, dphix] = eval_phi(x, basis);
   
    % Calculate control signal
    u = - 1/2 * inv(R) * gx' * dphix' * c_i + et;
    
else

    % *********************************************************************
    %
    % POST-LEARNING
    %
	% *********************************************************************    

    % Evaluate input gain matrix
    gx = eval_g(x, sys);    
    
    % Evaluate basis functions and gradient
    [~, dphix] = eval_phi(x, basis);    
    
    % Calculate control signal
    u = - 1/2 * inv(R) * gx' * dphix' * c_i;

end

