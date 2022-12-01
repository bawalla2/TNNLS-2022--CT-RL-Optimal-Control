function gx = eval_g(x, sys)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% EVALUATE SYSTEM INPUT GAIN MATRIX g(x)
%
% Brent Wallace  
%
% 2021-11-06
%
% This program, given a state vector x in R^n and system selection tag
% evaluates the input gain matrix g(x) corresponding to the example being
% simulated for.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% gx = eval_g(x, sys)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% x     Current value of the state (n-dimensional vector).
% sys   Tag for the system to be simulated for (string, see main.m).
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% gx    Input gain matrix g(x) evaluated at x corresponding to the system
%       sys (n x m matrix).
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
        
        gx =    [   0
                    sin(x(1))   ];
                
    % ***********************
    %
    % D. VRABIE, F.L. LEWIS (2009) -- HARD
    %
    
    case 'vrabie_lewis_2009_hard'
        
        gx =    [   0
                    sin(x(1))   ];
                                            
    % ***********************
    %
    % K.G. VAMVOUDAKIS, F.L. LEWIS (2010) -- F16 LINEAR
    %
    
    case 'vamvoudakis_lewis_2010_F16_lin'
        
        gx =    [   0
                    0
                    1   ];  

    % ***********************
    %
    % K.G. VAMVOUDAKIS, F.L. LEWIS (2010) -- 2ND ORDER NONLINEAR
    %     
    
    case 'vamvoudakis_lewis_2010_nonlin'
        
        gx =    [   0
                    cos(2 * x(1)) + 2 ];        

    % ***********************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE EXAMPLE
    %
    
    case 'jiang_jiang_2014_engine' 
        
        gx = -1;
                
    % ***********************
    %
    % T. BIAN, Z.-P. JIANG (2021) -- NONLINEAR EXAMPLE
    %
    % See Sec. V. A.
    %
    
    case 'bian_jiang_2021_nonlin_ex'
        
        % System parameters
        theta4 = 0.5; 
        
        gx = [  0
                theta4  ];
        

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
%     dxp = x(2);             % Cart velocity (m/s)
    theta = x(3);           % Pendulum angular displacement (rad)
%     dtheta = x(4);          % Pendulum angular velocity (rad/s)
    
    % g(x) -- \ddot{x}
    gxddx = 1 / (mc + mp * sin(theta)^2); 
    
    % g(x) -- \ddot{\theta}
    gxddtheta = - gxddx / l * cos(theta);

    % g(x)
    gx = [  0
            gxddx
            0
            gxddtheta   ];            
            
    % *********************************************************************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %   
    
    otherwise
        
        error('*** ERROR: SYSTEM TAG NOT RECOGNIZED ***');

end