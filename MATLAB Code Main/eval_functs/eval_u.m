function ux = eval_u(x, u_sett)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% EVALUATE FEEDBACK CONTROL LAW u(x)
%
% Brent Wallace  
%
% 2021-11-06
%
% This program, given a state vector x in R^n and preset tag evaluates the
% initial stabilizing policy u_0(x). 
%
% NOTE: If the initial policy can be implemented in terms of initial NN
% weights, then this function is not used. This function is used only when
% the initial policy varies from later iterate policies in a fashion which
% renders implementation impossible without treating the initial policy as
% a special case.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% ux = eval_u0(x, u_sett)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% x         Current value of the state (n-dimensional vector).
%
%           NOTE: In the case of an RADP unmatched algorithm, x will
%           consist of a (n + 1)-dimensional vector, partitioned as (cf.
%           eqns. (38)-(40)):
%   
%           x   : Known states (R^n)
%           z   : Known state where control acts (R^1)
%
%
% u_sett       (Struct) contains settings pertaining to the policy. Has the
%           following fields:
%   tag     (String) contains the string of the specific preset to evaluate
%           u_0(x) for.
%
% NOTE: If u_sett.tag = 'actor_g_known' (i.e., if the policy implemented is
% that of the actor network with g(x) known):
%
%           \hat{\mu}(x) = - 1/2 * R^{-1} g^T(x) \nabla \Phi^T(x) w
%
% then the following additional fields are needed in u_sett:
%   
%   basis   (Struct) basis structure
%   sys     (Struct) system used.
%   w       (N_1-dimensional vector) actor weights.
%   R       (m x m matrix) control penalty matrix.
% 
% If u_sett.tag = 'actor_g_known_lqr', then the following additional fields
% are needed in u_sett:
%
%   basis   (Struct) basis structure (contains weights, which have
%           been initialized programatically already).
%   sys     (Struct) system used.
%   R       (m x m matrix) control penalty matrix.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% ux    Initial stabilizing policy u_0(x), evaluated at x.
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

switch u_sett.tag

    % ***********************
    %
    % ZERO CONTROL SIGNAL
    %
    
    case '0'
        
        ux = 0;

    % ***********************
    %
    % ACTOR NETWORK WITH g(x) KNOWN
    %
    
    case 'actor_g_known'
        
        % Evaluate input dynamics g(x)
        gx = eval_g(x, u_sett.sys);
        
        % Evaluate basis functions
        [~, dphix] = eval_phi(x, u_sett.basis);
        
        % Evaluate policy
        ux = -1/2 * inv(u_sett.R) * gx' * dphix' * u_sett.w;

    % ***********************
    %
    % ACTOR NETWORK WITH g(x) KNOWN -- LQR
    %
    
    case 'actor_g_known_lqr'
        
        % Evaluate input dynamics g(x)
        gx = eval_g(x, u_sett.sys);
        
        % Evaluate basis functions
        [~, dphix] = eval_phi(x, u_sett.basis);
        
        % Evaluate policy
        ux = -1/2 * inv(u_sett.R) * gx' * dphix' * u_sett.basis.c0;        
        
    % ***********************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE EXAMPLE
    %
    % NOTE: For this example specifically, we have (cf. eqns.
    % (67)-(69))
    %
    % w = R
    % x = phi
    % z = psi
    %     
    
    case 'jiang_jiang_2014_engine_ex'
        
%         ux = 6 * x(1) - 2 * x(2);  % As shown in paper   
        ux = 9 * x(1) - 3 * x(2);   % As implemented in code
%     ux = 0;


    % ***********************
    %
    % T. BIAN, Z.-P. JIANG (2021) -- NONLINEAR EXAMPLE
    %

    case 'bian_jiang_2021_nonlin'

        k0 = 0;
        ux = -k0 * x(2);
        
    % ***********************
    %
    % COMPARISON D. VRABIE, F.L. LEWIS (2009) NONLINEAR "HARD" EXAMPLE
    % (SEC. 6.2)
    %     
    
    case 'comp_vrabie_lewis_2009_hard'
        
        ux = -0.5 * sin(x(1)) * ...
            (3 * x(2) - 0.2 * x(1)^2 * x(2) + 12 * x(2)^3);

    % ***********************
    %
    % D. VRABIE, F.L. LEWIS (2009) NONLINEAR "HARD" EXAMPLE (SEC. 6.2) --
    % STABILIZING POLICY IMPLEMENTABLE VIA THE MINIMUM CRITIC BASIS
    %     
    
    case 'vrabie_lewis_2009_hard_min_critic'
        
        ux = -0.5 * sin(x(1)) * ...
            (3 * x(2) + 12 * x(2)^3);        
        
    % ***********************
    %
    % D. VRABIE, F.L. LEWIS (2009) NONLINEAR "HARD" EXAMPLE (SEC. 6.2) --
    % STABILIZING POLICY IMPLEMENTABLE VIA sin(x_1) x_2 ONLY
    %     
    
    case 'vrabie_lewis_2009_hard_x_2_only'
        
        ux = -0.5 * sin(x(1)) * ...
            (10 * x(2)); 
        
    % ***********************
    %
    % D. VRABIE, F.L. LEWIS (2009) NONLINEAR "HARD" EXAMPLE (SEC. 6.2) --
    % TESTING
    %     
    
    case 'vrabie_lewis_2009_hard_testing'

        R = 1;
        
        gx = [ 0
                sin(x(1))   ];

        dphix = [   
                    2 * x(1)            0
                    0                   2 * x(2)
%                     0                   4 * x(2)^3  % Opt. functs. end here  
                    x(2)                x(1)  
                    4 * x(1)^3          0
                    3 * x(1)^2 * x(2)   x(1)^3
                    2 * x(1) * x(2)^2   2 * x(1)^2 * x(2)
                    x(2)^3              3 * x(1) * x(2)^2                    
                                                            ];          
                                
        c = [   0.1871
                1.8123
                0
                0.2025
                0
                0.1876
                0       ];
    
        ux = - 1/2 * inv(R) * gx' * dphix' * c;
        
    % ***********************
    %
    % COMPARISON CART INVERTED PENDULUM
    %     
    
    case 'comp_cip'
        
        % State feedback u = - K x obtained from LQR
        ux = - [ -31.6228  -50.2741 -247.4658  -63.6293] * x;     

    
    % *********************************************************************
    % *********************************************************************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %  
    % *********************************************************************
    % *********************************************************************
    
    otherwise
        
        error('*** ERROR: PRESET TAG NOT RECOGNIZED ***');

end