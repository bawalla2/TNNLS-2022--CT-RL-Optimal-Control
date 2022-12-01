function [sys, sys_plot_settings] = config_sys(sys)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INITIALIZE SYSTEM PARAMETERS
%
% Brent Wallace  
%
% 2022-02-08
%
% This program, given a system tag, initializes system parameters such as
% system order, number of inputs, etc. See config_main.m for a list of
% system options.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% sys = config_sys(sys)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% sys           (Struct) System for design. Has the following fields:
%   tag         (String) System tag (see above for options).
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% sys           (Struct) Contains system parameters, fully initialized. Has
%               the following fields:
%   tag         (String) System tag (see above for options).
%   n           (Integer) System order.
%   m           (Integer) Number of inputs of system.
%   total_order (Integer) Total system order (in most cases, this is just
%               n, but for e.g. RADP with matched uncertainty it may be n
%               + p, where p denotes the order of the unknown dynamics).
%
% sys_plot_settings    (Struct) Contains plot settings corresponding to
%                       the system designed for. Has the following fields:
%   x_sclvec    (n-dim. Vector) A vector containing the scaling desired for
%               each state variable on plots. Is a vector of ones by
%               default (i.e., no scaling). E.g., if x_2(t) is an angle in
%               radians, but it is desired to plot it in degrees, declare
%               x_sclvec(2) = 180/pi.
%   x_t_title   (n-dim. Cell) i-th entry contains the title to put
%                       on the plot of x_i(t). E.g., x_t_title{i} =
%                       'Pendulum angle \theta'.
%   x_t_xlabel  (n-dim. Cell) i-th entry contains the x-axis label
%                       to put on the plot of x_i(t). E.g., x_t_xlabel{i} =
%                       '\theta(t) (deg)'.
%   x_t_filename (n-dim. Cell) i-th entry contains the desired file
%                       name for the plot of x_i(t). E.g., x_t_filename{2}
%                       = 'x_2_t' or 'theta_t'.
%                       to put on the plot of x_i(t). E.g., x_t_xlabel{i} =
%                       '\theta(t) (deg)'.
%   u_t_title   (m-dim. Cell) i-th entry contains the title to put
%                       on the plot of u_i(t). E.g., u_t_title{i} = 'Thrust
%                       T'. For m = 1, this could simply be, e.g.,
%                       u_t_title{1} = 'Control Signal u(t)'
%   u_t_xlabel  (m-dim. Cell) i-th entry contains the x-axis label
%                       to put on the plot of u_i(t). E.g., u_t_xlabel{i} =
%                       'T(t)'. For m = 1, this could simply be, e.g.,
%                       u_t_xlabel{1} = 'u(t)'.
%           u_t_filename (m-dim. Cell) i-th entry contains the desired file
%                       name for the plot of u_i(t). E.g., u_t_filename{2}
%                       = 'u_2_t' or 'T_t'.
%
% *************************************************************************
% *************************************************************************
% *************************************************************************



%%
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
% *************************************************************************
%
% CONFIGURE SYSTEM PARAMETERS
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

switch sys.tag

    % ***********************
    %
    % D. VRABIE, F.L. LEWIS (2009) -- EASY
    %
    
    case 'vrabie_lewis_2009_easy'
        
        % System parameters
        sys.n = 2;                  % System order
        sys.m = 1;                  % Number of inputs   
        sys.total_order = sys.n;    % Total system order
        
    % ***********************
    %
    % D. VRABIE, F.L. LEWIS (2009) -- HARD
    %
    
    case 'vrabie_lewis_2009_hard'
        
        % System parameters
        sys.n = 2;                  % System order
        sys.m = 1;                  % Number of inputs  
        sys.total_order = sys.n;    % Total system order
        
    % ***********************
    %
    % K.G. VAMVOUDAKIS, F.L. LEWIS (2010) -- F16 LINEAR
    %
    
    case 'vamvoudakis_lewis_2010_F16_lin'
        
        % System parameters
        sys.n = 3;                  % System order
        sys.m = 1;                  % Number of inputs
        sys.total_order = sys.n;    % Total system order
        
           
    % ***********************
    %
    % K.G. VAMVOUDAKIS, F.L. LEWIS (2010) -- 2ND ORDER NONLINEAR
    %
    
    case 'vamvoudakis_lewis_2010_nonlin'
        
        % System parameters
        sys.n = 2;                  % System order
        sys.m = 1;                  % Number of inputs    
        sys.total_order = sys.n;    % Total system order
        
    % ***********************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE
    %
    
    case 'jiang_jiang_2014_engine'
        
        % System parameters
        sys.p = 1;                          % Order of unknown w-dynamics
        sys.n = 1;                          % Order of known x-dynamics
        sys.m = 1;                          % Number of inputs    
        sys.total_order = sys.n + sys.p + 1;    % Total system order
        
    % ***********************
    %
    % T. BIAN, Z.-P. JIANG (2021) -- NONLINEAR EXAMPLE
    %
    
    case 'bian_jiang_2021_nonlin_ex'
        
        % System parameters
        sys.n = 2;                  % Order of system   
        sys.m = 1;                  % Number of inputs
        sys.total_order = sys.n;    % Total system order
 
    % ***********************
    %
    % CART INVERTED PENDULUM
    %
    
    case 'cip'
        
        % System parameters
        sys.n = 4;                  % Order of system   
        sys.m = 1;                  % Number of inputs      
        sys.total_order = sys.n;    % Total system order
        
        % System constants
        mc = 1;         % Cart mass (kg)
        mp = 0.1;       % Pendulum ball mass (kg)
        l = 0.5;        % Pendulum length (m)
        g = 9.81;       % Gravitational field constant (m/s^2)
        
        % Linearization "A" matrix
        A =  [  0       1       0               0
                0       0       -(mp*g/mc)        0
                0       0       0               1
                0       0       (g/l)*(1+(mp/mc)) 0  ];
        
        
        % Store system constants
        sys.mc = mc;
        sys.mp = mp;
        sys.l = l;
        sys.g = g;
        
        % Store linearization "A" matrix
        sys.A = A;
        
        
    % ***********************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %   
    
    otherwise
        
        error('*** ERROR: SYSTEM TAG NOT RECOGNIZED ***');

end


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CONFIGURE SYSTEM PLOT SETTINGS
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% *************************************************************************
%
% DEFAULT SETTINGS
%
% Note: These settings will be initialized as default. If different
% settings are desired, simply overwrite the default value in the switch
% case structure below
%
% *************************************************************************

% Initialize empty structs
x_t_title = cell(sys.n, 1);
x_t_ylabel = cell(sys.n, 1);
x_t_filename = cell(sys.n, 1);
x_t_state_lgd = cell(sys.n, 1);
u_t_title = cell(sys.m, 1);
u_t_ylabel = cell(sys.m, 1); 
u_t_filename = cell(sys.m, 1);

% State variable scaling for plots (default, no scaling)
x_sclvec = ones(sys.n, 1);

% Fill out state trajectory settings
if sys.n == 1
   
    % n = 1. Do not include subscripts in plots.
    x_t_title{1} = ['State Trajectory $x(t)$'];
    x_t_ylabel{1} = ['$x(t)$'];
    x_t_state_lgd{1} = x_t_ylabel{1};
    x_t_filename{1} = ['x_t'];
    
else
    
    % n > 1. Include subscripts in plots.
    for i = 1:sys.n  
        nsi = num2str(i);
        x_t_title{i} = ['State Trajectory $x_{' nsi '}(t)$'];
        x_t_ylabel{i} = ['$x_{' nsi '}(t)$'];  
        x_t_state_lgd{i} = x_t_ylabel{i};
        x_t_filename{i} = ['x_' nsi '_t'];
    end
    
end

% Fill out control signal settings
if sys.m == 1
   
    % m = 1. Do not include subscripts in plots.
    u_t_title{1} = 'Control Signal $u(t)$';
    u_t_ylabel{1} = '$u(t)$';  
    u_t_filename{1} = 'u_t';
    
else
    
    % m > 1. Include subscripts in plots.
    for i = 1:sys.m   
        nsi = num2str(i);
        x_t_title{i} = ['Control Signal $u_{' nsi '}(t)$'];
        x_t_ylabel{i} = ['$u_{' nsi '}(t)$'];  
        x_t_filename{i} = ['u_' nsi '_t'];
    end    
    
end


% *************************************************************************
%
% CUSTOM SETTINGS
%
% Note: If default settings are desired, simply leave the respective system
% case empty.
%
% *************************************************************************

switch sys.tag

    % ***********************
    %
    % D. VRABIE, F.L. LEWIS (2009) -- EASY
    %
    
    case 'vrabie_lewis_2009_easy'

        % Use default settings     
        
    % ***********************
    %
    % D. VRABIE, F.L. LEWIS (2009) -- HARD
    %
    
    case 'vrabie_lewis_2009_hard'
        
        % Use default settings       
        
    % ***********************
    %
    % K.G. VAMVOUDAKIS, F.L. LEWIS (2010) -- F16 LINEAR
    %
    
    case 'vamvoudakis_lewis_2010_F16_lin'
        
        % Use default settings
           
    % ***********************
    %
    % K.G. VAMVOUDAKIS, F.L. LEWIS (2010) -- 2ND ORDER NONLINEAR
    %
    
    case 'vamvoudakis_lewis_2010_nonlin'
        
        % Use default settings      
        
    % ***********************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE
    %
    
    case 'jiang_jiang_2014_engine'

        % Overwrite state trajectory plot settings
        x_t_title = cell(3, 1);
        x_t_ylabel = cell(3, 1);
        x_t_filename = cell(3, 1);
        
        % Fill values
        x_t_title{1} = 'Rotating Stall Amplitude R(t)';
        x_t_ylabel{1} = 'R(t)';
        x_t_filename{1} = 'R_t';
        x_t_title{2} = ['Scaled Annulus-Averaged Flow \phi(t)' ...
                        ' = \Phi(t) - \Phi_e'];
        x_t_ylabel{2} = '\phi(t)';
        x_t_filename{2} = 'phi_t';
        x_t_title{3} = ['Plenum Pressure Rise \psi(t)' ...
                        ' = \Psi(t) - \Psi_e'];
        x_t_ylabel{3} = '\psi(t)';
        x_t_filename{3} = 'psi_t';
        
        
    % ***********************
    %
    % T. BIAN, Z.-P. JIANG (2021) -- NONLINEAR EXAMPLE
    %
    
    case 'bian_jiang_2021_nonlin_ex'
        
        % Use default settings 

    % ***********************
    %
    % CART INVERTED PENDULUM
    %
    
    case 'cip'

        % Change scaling in x_3 = \theta, x_4 = \dot{\theta} from radians
        % to degrees
%         x_sclvec(3:4) = 180/pi;
        
        % Overwrite state trajectory plot settings
        x_t_title = cell(4, 1);
        x_t_ylabel = cell(4, 1);
        x_t_filename = cell(4, 1);
        
        % Fill values
        x_t_title{1} = 'Cart Position x(t)';
        x_t_ylabel{1} = '$x(t)$ (m)';
        x_t_state_lgd{1} = '$x(t)$';
        x_t_filename{1} = 'x_t';
        
        x_t_title{2} = 'Cart Velocity $\dot{x}$(t)';
        x_t_ylabel{2} = '$\dot{x}(t)$ (m/s)';
        x_t_state_lgd{2} = '$\dot{x}(t)$';
        x_t_filename{2} = 'dx_t';
        
        x_t_title{3} = 'Pendulum Angular Displacement $\theta$(t)';
%         x_t_ylabel{3} = '$\theta$(t) (deg)';
        x_t_ylabel{3} = '$\theta(t)$ (rad)';
        x_t_state_lgd{3} = '$\theta(t)$';
        x_t_filename{3} = 'theta_t'; 
        
        x_t_title{4} = 'Pendulum Angular Velocity $\dot{\theta}$(t)';
        x_t_ylabel{4} = '$\dot{\theta}(t)$ (deg/s)';
        x_t_state_lgd{4} = '$\dot{\theta}(t)$';
        x_t_filename{4} = 'dtheta_t'; 
        
        
    % ***********************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %   
    
    otherwise
        
        error('*** ERROR: SYSTEM TAG NOT RECOGNIZED ***');

end

% *************************************************************************
%
% STORE SYSTEM PLOT SETTINGS
% 
% *************************************************************************

sys_plot_settings.x_sclvec = x_sclvec;

sys_plot_settings.x_t_title = x_t_title;
sys_plot_settings.x_t_ylabel = x_t_ylabel;
sys_plot_settings.x_t_state_lgd = x_t_state_lgd;
sys_plot_settings.x_t_filename = x_t_filename;
sys_plot_settings.u_t_title = u_t_title;
sys_plot_settings.u_t_ylabel = u_t_ylabel;
sys_plot_settings.u_t_filename = u_t_filename;

% *************************************************************************
%
% STORE SYSTEM LINEARIZATION 'B" MATRIX
% 
% *************************************************************************

sys.B = eval_g(zeros(sys.n,1), sys);
