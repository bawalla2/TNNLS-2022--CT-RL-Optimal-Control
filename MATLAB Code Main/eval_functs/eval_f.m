function fx = eval_f(x, sys)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% EVALUATE SYSTEM DRIFT DYNAMICS
%
% Brent Wallace  
%
% 2021-11-06
%
% This program, given a state vector x in R^n and system selection tag
% evaluates the drift dynamics f(x) corresponding to the example being
% simulated for.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% fx = eval_f(x, sys)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% x     Current value of the state (n-dimensional vector).
%
%       NOTE: In the case of an RADP matched algorithm (when the system
%       dynamic uncertainty is included), x will consist of a (p +
%       n)-dimensional vector, partitioned as (cf. eqns. (8)-(9)):
%   
%       w   : Uncertain states (R^p)
%       x   : Known states (R^n)
%
%       NOTE: In the case of an RADP unmatched algorithm, x will consist of
%       a (p + n + 1)-dimensional vector, partitioned as (cf. eqns.
%       (38)-(40)):
%   
%       w   : Uncertain states (R^p)
%       x   : Known states (R^n)
%       z   : Known state where control acts (R^1)
%
%
% sys   Array consisting of the following fields:
%
%   tag         (String) System tag.
%   n           (Integer) System order n.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% fx    Drift dynamics f(x) evaluated at x corresponding to the system sys
%       (n-dimensional vector).
%
%       NOTE: In the case of an RADP unmatched algorithm, fx will consist
%       of a (p + n + 1 + 1 + 1)-dimensional vector partitioned as (cf.
%       eqns. (38)-(40)):
%
%           Delta_w(w, x)   : w-dynamic uncertainty (R^p)
%           f(x)            : Nominal x-drift dynamics (R^n)
%           Delta(w, x)     : x-dynamic uncertainty (R^1)
%           f_1(x, z)       : Nominal z-drift dynamics (R^1)
%           Delta_1(w, x, z): z-dynamic uncertainty (R^1)
%
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% BEGIN MAIN
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
        
        fx =    [     -1     1
                    -1/2    -1/2 ] * x + ...
                [   0
                    1/2*x(2)*(sin(x(1)))^2   ];
                
    % ***********************
    %
    % D. VRABIE, F.L. LEWIS (2009) -- HARD
    %
    
    case 'vrabie_lewis_2009_hard'
        
        fx =    [     -1     1
                    -1/2    -1/2 ] * x + ...
                [   2*(x(2))^3
                    1/2 * x(2) * (1 + 2 * x(2)^2) * (sin(x(1)))^2   ];
                                            
    % ***********************
    %
    % K.G. VAMVOUDAKIS, F.L. LEWIS (2010) -- F16 LINEAR
    %
    
    case 'vamvoudakis_lewis_2010_F16_lin'
        
        fx =    [   -1.01887    0.90506     -0.00215
                    0.82225     -1.07741    -0.17555
                    0           0           -1          ] * x;  

    % ***********************
    %
    % K.G. VAMVOUDAKIS, F.L. LEWIS (2010) -- 2ND ORDER NONLINEAR
    %     
    
    case 'vamvoudakis_lewis_2010_nonlin'

        fx =    [     -1     1
                    -1/2    -1/2 ] * x + ...
                [   0
                    1/2 * x(2) * (cos(2 * x(1)) + 2)^2   ];        
        
    % ***********************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE EXAMPLE
    %
    % NOTE;
    %
    % For RADP unmatched examples, the input argument "x" and output
    % argument "fx" are partitioned as described above
    %
    % Note that for this example specifically, we have (cf. eqns.
    % (67)-(69))
    %
    % w = R
    % x = phi
    % z = psi
    %     
    
    case 'jiang_jiang_2014_engine'
        
        % System parameters (cf. Sec. V. A.)
        sigma = 0.3;    
        beta = 0.702;
        Psi_C0 = 1.7;
        
        % State dimensions
        n = sys.n;              % Order of known x-dynamics
        p = sys.p;       % Order of unknown w-dymanics
        
        % Extract states
        w = x(1:p);             % Unknown w-dynamics
        xs = x(p+1:p+n);        % Known x-dynamics
        z = x(p+n+1);           % Known z-dynamics

        % Evaluate Psi_C(phi + 1)
        Psi_C_xp1 = Psi_C0 + 2 - 1.5 * xs^2 - 0.5 * xs^3;

        % Evaluate dynamics (cf. eqns. (67)-(69))
        %
        %           Delta_w(w, x)   : w-dynamic uncertainty (R^p)
        %           f(x)            : Nominal x-drift dynamics (R^n)
        %           Delta(w, x)     : x-dynamic uncertainty (R^n)
        %           f_1(x, z)       : Nominal z-drift dynamics (R^1)
        %           Delta_1(w, x, z): z-dynamic uncertainty (R^1)
        
        Delta_w = - sigma * w^2 - sigma * w * (2 * xs + xs^2);
        f_x = Psi_C_xp1 - 2 - Psi_C0;
        Delta = 3 * w * xs + 3 * w;
        f1_xz = 0;
        Delta_1 = 0;
        
        % Pack dynamics into output argument "fx"       
        fx =    [   Delta_w
                    f_x
                    Delta
                    f1_xz
                    Delta_1     ];
        
    % ***********************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- POWER SYSTEM EXAMPLE
    %
    % NOTE;
    %
    % For RADP unmatched examples, the input argument "x" and output
    % argument "fx" are partitioned as described above
    %
    % Note that for this example specifically, we have (cf. eqns.
    % (79)-(82))
    %
    % w
    % x1
    % x_2
    % z
    %     
    
    case 'jiang_jiang_2014_power_sys'

        % System parameters (cf. Sec. V. B.)
        D = 5;
        H = 4;              % Inertia constant
        omega_0 = 314.159;  
        x_T = 0.127;        % Transformer reactance
        x_L = 0.4853;       % Line reactance
        x_d = 1.863;        % Direct axis reactance
        x_dp = 0.257;       % Direct axis transient reactance
        T_d0p = 0.5;        % Direct axis transient time constant
        delta_0 = 1.2566;   % Steady state rotor angle
        Vs = 1;             % Voltage on infinite bus
        T_T = 2;            
        K_G = 1;
        K_T = 1;
        T_G = 0.2;
        
        % Constants dependent on above parameters (cf. Sec. V. B.)
        x_ds = x_T + x_L + x_d;
        x_dsp = x_T + x_L + x_dp;
        
        a_1 = 1 / T_d0p;                % ???
        a_2 = (x_d - x_dp) / T_d0p * Vs^2 / (x_dsp * x_ds); % ???
        a_3 = delta_0;
        
        b_1 = D / (2 * H);
        b_2 = (omega_0 / H) * (Vs / x_ds);
        
        
        % State dimensions
        n = sys.n;              % Order of known x-dynamics
        p = sys.p;       % Order of unknown w-dymanics
        
        % Extract states
        w = x(1:p);             % Unknown w-dynamics
        xs = x(p+1:p+n);        % Known x-dynamics
        z = x(p+n+1);           % Known z-dynamics

        % Evaluate Psi_C(phi + 1)
        Psi_C_xp1 = Psi_C0 + 2 - 1.5 * xs^2 - 0.5 * xs^3;

        % Evaluate dynamics (cf. eqns. (67)-(69))
        %
        %           Delta_w(w, x)   : w-dynamic uncertainty (R^p)
        %           f(x)            : Nominal x-drift dynamics (R^n)
        %           Delta(w, x)     : x-dynamic uncertainty (R^n)
        %           f_1(x, z)       : Nominal z-drift dynamics (R^1)
        %           Delta_1(w, x, z): z-dynamic uncertainty (R^1)
        
        Delta_w = - sigma * w^2 - sigma * w * (2 * xs + xs^2);
        f_x = Psi_C_xp1 - 2 - Psi_C0;
        Delta = 3 * w * xs + 3 * w;
        f1_xz = 0;
        Delta_1 = 0;
        
        % Pack dynamics into output argument "fx"       
        fx =    [   Delta_w
                    f_x
                    Delta
                    f1_xz
                    Delta_1     ];
            
    
    % ***********************
    %
    % T. BIAN, Z.-P. JIANG (2021) -- NONLINEAR EXAMPLE
    %
    % See Sec. V. A.
    %
    
    case 'bian_jiang_2021_nonlin_ex'
        
        % System parameters
        theta1 = -0.5;
        theta2 = -2;
        theta3 = 0.125;
        %theta4 = 0.5;  % Not needed here
        
        fx = [  theta1 * x(1)^3 - x(1) + theta2 * x(2)
                theta3 * x(2)^3 - x(2)                  ];
                
 
    % ***********************
    %
    % CART INVERTED PENDULUM
    %
    
    case 'cip'
    
    % System constants
    mc = sys.mc;    % Cart mass (kg)
    mp = sys.mp;    % Pendulum ball mass (kg)
    l = sys.l;      % Pendulum length (m)
    g = sys.g;      % Gravitational field constant (m/s^2)
    
    % Unpack state vector
%     xp = x(1);              % Cart position (m)
    dxp = x(2);             % Cart velocity (m/s)
    theta = x(3);           % Pendulum angular displacement (rad)
    dtheta = x(4);          % Pendulum angular velocity (rad/s)
    
    % x dynamics
    ddxp = (mp * l * dtheta^2 * sin(theta) ...
                - mp * g * sin(theta) * cos(theta)) ...
                / (mc + mp * sin(theta)^2);
            
    % theta dynamics
    ddtheta = - ddxp / l * cos(theta) + g / l * sin(theta);
    
    % Pack output vector
    fx = [  dxp
            ddxp
            dtheta
            ddtheta     ];
    
               
            
                
    % ***********************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %   
    
    otherwise
        
        error('*** ERROR: SYSTEM TAG NOT RECOGNIZED ***');

end