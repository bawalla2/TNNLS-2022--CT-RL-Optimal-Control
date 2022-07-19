function [phix, dphix] = eval_phi(x, basis)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% EVALUATE BASIS FUNCTIONS AND BASIS FUNCTION GRADIENTS
%
% Brent Wallace  
%
% 2021-11-06
%
% This program, given a state vector x in R^n and N-dimensional basis
% selection tag evaluates the basis functions and their gradients
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% [phix, dphix] = eval_phi(x, basis)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% x         Current value of the state (n-dimensional vector).
% basis     (Struct) Contains basis parameters. Includes the fields:
%   tag     (String) Tag of the corresponding basis.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% phix      Basis functions evaluated at x (N-dimensional vector)
% dphix     Basis fuction gradient evaluated at x (N x n matrix)
%
% *************************************************************************
%
% NOTE ON 'monomial' BASIS
%
% *************************************************************************
%
% If the basis is of tag 'monomial', the following additional fields must
% be present in the 'basis' struct:
%
%   type        (String) Monomial basis type. Has the following options:
%       'all'       Evaluate all monomials of total degree <= K.
%       'even'      Evaluate even monomials of total degree <= K.
%       'odd'       Evaluate odd monomials of total degree <= K.
%       'total_deg_K' Monomials of total degree = K only.
%       'custom'    Evaluate custom list of monomials.
%
%   IC_type     (String) Describes the network weight initialization
%               behavior. Has the following options:
%       'none'      ICs are to be declared manually later. Do not
%                   initialize weights here.
%       'custom'    ICs are to be set programatically here. See below note.
%       'lqr'       This network is to be implemented in an actor network
%                   with g(x) known:
%           \hat{\mu}(x) = - 1/2 R^{-1} g^T(x) \nabla \Phi^T(x) w
%                   and the weights w_0 \in R^{N_1} are to be initialized
%                   according to the procedure outlined in
%                   config_basis_ICs_lqr.m. See this function description
%                   for more details.
%
% If 'type' is 'all', 'even', or 'odd', the following additional field
% needs to be declared in the 'basis' struct:
%   
%   K           (Integer) Max total degree of monomials in basis.
%
% If 'type' is 'custom', the following additional field needs to be
% declared in the 'basis' struct:
%
%   dmat        (N x n Matrix) A matrix of integer elements, where N
%               denotes the total number of activation functions and n the
%               system order. The i-th row, j-th column of of dmat will
%               contain the integer corresponding to the power of x_j in
%               the activation function \phi_j(x). E.g., for n = 3, j = 2,
%               dmat(j,:) might look like [ 0 2 1 ]. Then \phi_j(x) = x_1^0
%               * x_2^2 * x_3^1.
%
% NOTE: If 'type' is 'all', 'even', or 'odd', the field 'dmat' is
% initialized automatically in config.m.
%
% NOTE: If 'IC_type' is 'custom', then the function config_basis_ICs.m will
% be executed to initialize the network weights. The following additional
% fields need to be declared 'basis' struct (see config_basis_ICs.m for a
% description of these fields):
%   
%   dmat_IC     ('numIC' x n Matrix)
%   c0_ic       ('numIC'-dimensional vector)
%   
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

switch basis.tag

    
    % *********************************************************************
    %
    % MONOMIAL BASIS
    %
    % Note: See config.m for details of how the basis was initialized.
    % Depending on the value in 'basis.type', this will evaluate the
    % following monomials:
    %
    %   'all'       All monomials of total degree <= K.
    %   'even'      Even monomials of total degree <= K.
    %   'odd'       Odd monomials of total degree <= K.
    %   
    
    case 'monomial'
        
        % ***********************
        %
        % EXTRACT BASIS PARAMETERS
        %
        
        % Holds monomial powers of activation functions
        dmat = basis.dmat;
        
        % Number of activation functions
        N = size(dmat, 1);
        
        % System order
        n = size(dmat, 2);   
        
        % ***********************
        %
        % Phi(x) CALCULATION
        %
        % [Phi(x)](i) = \phi_i(x)
        %
        %             = \prod_{k=1}^{n} x_k^{dmat(i,k)}
        %
        
        phix = ones(N, 1);
        
        for i = 1:N
           
            for k = 1:n
               
                phix(i) = phix(i) * x(k)^dmat(i,k);
                
            end
            
        end
        
        % ***********************
        %
        % \nabla Phi(x) CALCULATION
        %
        % [\nabla Phi(x)](i,j) = \partial \phi_i(x) / \partial x_j
        %
        %             = [ \prod_{k \neq j} x_k^{dmat(i,k)} ] * 
        %                  dmat(i,j) * x_j^{dmat(i,k) - 1}
        %
        
        dphix = ones(N, n);
        
%         for i = 1:N
%             for j = 1:n
%                 % Calculate [\nabla Phi(x)](i,j)
%                 for k = 1:n
%                     % Handle x_k = 0 case
%                     if abs(x(k)) == 0
%                         dphix(i,j) = 0;
%                     % Else, x_k ~= 0
%                     else    
%                         if k == j 
%                             dphix(i,j) = dphix(i,j) * ...
%                                 (dmat(i,k)) * x(k)^(dmat(i,k) - 1);  
%                         else
%                             dphix(i,j) = dphix(i,j) * x(k)^(dmat(i,k));
%                         end
%                     end
%                 end 
%             end
%         end
        
        
        for i = 1:N
            for j = 1:n
                % Calculate [\nabla Phi(x)](i,j)
                % Handle x_j = 0 case
                if abs(x(j)) == 0
                    dphix(i,j) = 0;
                % Else, x_j ~= 0
                else   
                    dphix(i,j) = (dmat(i,j)) * phix(i) / x(j);
                end
            end
        end
        
    
    % *********************************************************************
    %
    % 1ST ORDER SYSTEM, DEGREE 4
    %
    
    case 'order_1_degree_4'
        
        phix = [    x(1)
                    x(1)^2
                    x(1)^3
                    x(1)^4        ];
        
        dphix = [   1      
                    2*x(1)     
                    3*x(1)^2
                    4*x(1)^3        ];    
    
    % *********************************************************************
    %
    % 2ND ORDER SYSTEM, DEGREE 2
    %
    
    case 'order_2_degree_2'
        
        phix = [    x(1)^2
                    x(1) * x(2)
                    x(2)^2        ];
        
        dphix = [   2 * x(1)    0      
                    x(2)        x(1)     
                    0           2 * x(2)    ];       
                
    % *********************************************************************
    %
    % 2ND ORDER SYSTEM, DEGREE 4
    %
    
    case 'order_2_degree_4'
        
        phix = [    x(1)^2
                    x(1) * x(2)
                    x(2)^2
                    x(1)^4
                    x(1)^3 * x(2)
                    x(1)^2 * x(2)^2
                    x(1) * x(2)^3
                    x(2)^4              ];
        
        dphix = [   2 * x(1)            0
                    x(2)                x(1)
                    0                   2 * x(2)
                    4 * x(1)^3          0
                    3 * x(1)^2 * x(2)   x(1)^3
                    2 * x(1) * x(2)^2   2 * x(1)^2 * x(2)
                    x(2)^3              3 * x(1) * x(2)^2
                    0                   4 * x(2)^3          ];               
                
    % *********************************************************************
    %
    % 2ND ORDER SYSTEM, DEGREE 4 -- ALL TERMS
    %
    % Let i index x(1) (1 <= i <= 5), j index x(2) (1 <= j <= 5)
    % 
    % phix(5(i-1)+(j-1)+1) = x(1)^(i-1) * x(2)^(j-1)
    %
    % dphix(5(i-1)+(j-1)+1,:) = 
    %   [   (i-1)x(1)^(i-2)*x(2)^(j-1)  x(1)^(i-1)*(j-1)x(2)^(j-2)  ]
    
    case 'order_2_degree_4_all_terms'
        
        phix = zeros(25,1);
        dphix = zeros(25,2);
        
        for i = 1:5
            
            for j = 1:5
                
                phix(5*(i-1)+(j-1)+1) = x(1)^(i-1) * x(2)^(j-1);
                
                dphix(5*(i-1)+(j-1)+1,1) = (i-1)*x(1)^(i-2) * x(2)^(j-1);  
                
                dphix(5*(i-1)+(j-1)+1,2) = x(1)^(i-1) * (j-1)*x(2)^(j-2);
                
            end
                     
        end
        
        % Remove zeroth-order terms
        phix = phix(2:end);
        dphix = dphix(2:end,:);
                                            
    % *********************************************************************
    %
    % 3RD ORDER SYSTEM, DEGREE 2
    %
    
    case 'order_3_degree_2'
        
        phix = [    x(1)^2
                    x(1) * x(2)
                    x(1) * x(3)
                    x(2)^2
                    x(2) * x(3)
                    x(3)^2          ];
        
        dphix = [   2 * x(1)    0           0
                    x(2)        x(1)        0
                    x(3)        0           x(1)
                    0           2 * x(2)    0
                    0           x(3)        x(2)
                    0           0           2 * x(3)    ];
 

    % *********************************************************************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE EXAMPLE -- Phi1
    %
    % Critic NN activation function parameters (cf. eqn. (13))
    %       
    
    case 'jiang_jiang_2014_engine_ex_Phi1'
    
        phix = [    x^2
                    x^3
                    x^4     ];
                
        dphix = [];     % Unused
    
    % *********************************************************************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE EXAMPLE -- Phi2
    %
    % Actor NN activation function parameters (cf. eqn. (14))
    % 
    
    case 'jiang_jiang_2014_engine_ex_Phi2'
    
        phix = [    x
                    x^2
                    x^3     ];                   
        
        dphix = [];     % Unused
        
    % *********************************************************************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE EXAMPLE -- Phi3
    %
    % NN activation function parameters for approximation of \hat{f}_1(x,z)
    % (cf. eqn. (43))
    %     
    
    case 'jiang_jiang_2014_engine_ex_Phi3'
        
        phix = [    x(1)
                    x(1)^2
                    x(1)^3
                    x(2)
                    x(1) * x(2)
                    x(1)^2 * x(2)
                    x(1)^3 * x(2)   ];
                           
        dphix = [];     % Unused
 
    % *********************************************************************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE EXAMPLE -- Phi4
    %
    % NN activation function parameters for approximation of \hat{g}_1(x)
    % (cf. eqn. (44))
    %     
    
    case 'jiang_jiang_2014_engine_ex_Phi4'
        
        phix = [    1
                    x
                    x^2
                    x^3     ];         
        
        dphix = [];     % Unused

        
    % *********************************************************************
    %
    % T. BIAN, Z.-P. JIANG (2021) -- NONLINEAR EXAMPLE HAMILTONIAN NN
    % \Theta(x)
    %        
    
    case 'bian_jiang_2021_nonlin_theta'
        
        x1 = x(1);
        x2 = x(2);


        phix = [    x1
                    x2
                    x1^2
                    x1 * x2
                    x2^2
                    x1^3
                    x1^2 * x2
                    x1 * x2^2
                    x2^3
                    x1^4
                    x1^3 * x2
                    x1^2 * x2^2
                    x1 * x2^3
                    x2^4          ];
        
                
        dphix = [];     % Unused     

    % *********************************************************************
    %
    % T. BIAN, Z.-P. JIANG (2021) -- NONLINEAR EXAMPLE HAMILTONIAN NN
    % \Psi_1(x)
    %        
    
    case 'bian_jiang_2021_nonlin_psi'

        phix = [    x(1)
                    x(2)  ];
        
                
        dphix = [];     % Unused             
         
         
    % *********************************************************************
    %
    % ALL METHOD COMPARISON OF D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- BASIS FOR CRITIC WITH MINIMAL NUMBER OF
    % ACTIVATION FUNCTIONS
    %             
    
    case 'comp_vrabie_lewis_2009_hard_critic_min'  
        
        
        phix = [    x(1)^2
                    x(2)^2
                    x(2)^4              ];
        
        dphix = [   2 * x(1)            0
                    0                   2 * x(2)
                    0                   4 * x(2)^3          ];           

    % *********************************************************************
    %
    % D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- MINIMAL BASIS FOR CRITIC, WITH SOME
    % ACTIVATION FUNCTIONS ADDED -- N_1 = 4 TERMS
    %             
    
    case 'vrabie_lewis_2009_hard_critic_nonmin'  
        
        
        phix = [    x(1)^2
                    x(2)^2
                    x(2)^4              
                    x(1) * x(2)     ];
        
        dphix = [   2 * x(1)            0
                    0                   2 * x(2)
                    0                   4 * x(2)^3
                    x(2)                x(1)            ];    
                
                
    % *********************************************************************
    %
    % D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- MINIMAL BASIS FOR CRITIC, WITH SOME
    % ACTIVATION FUNCTIONS ADDED -- N_1 = 8 TERMS
    %             
    
    case 'vrabie_lewis_2009_hard_critic_8terms'  
        
        
        phix = [    
                    x(1)^2
                    x(2)^2
                    x(2)^4                  % Opt. functs. end here  
                    x(1) * x(2) 
                    x(1)^4
                    x(1)^3 * x(2)
                    x(1)^2 * x(2)^2
                    x(1) * x(2)^3                    
                                    ];
        
        dphix = [   
                    2 * x(1)            0
                    0                   2 * x(2)
                    0                   4 * x(2)^3  % Opt. functs. end here  
                    x(2)                x(1)  
                    4 * x(1)^3          0
                    3 * x(1)^2 * x(2)   x(1)^3
                    2 * x(1) * x(2)^2   2 * x(1)^2 * x(2)
                    x(2)^3              3 * x(1) * x(2)^2                    
                                                            ];     
                                                        
                                                        
    % *********************************************************************
    %
    % D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- MINIMAL BASIS FOR CRITIC, WITH SOME
    % ACTIVATION FUNCTIONS ADDED -- N_1 = 7 TERMS, ESSENTIAL ACTIVATION
    % FUNCTION x_2^4 ABSENT
    %             
    
    case 'vrabie_lewis_2009_hard_critic_7terms_no_x_2^4'  
        
        
        phix = [    
                    x(1)^2
                    x(2)^2
%                     x(2)^4                  % Opt. functs. end here  
                    x(1) * x(2) 
                    x(1)^4
                    x(1)^3 * x(2)
                    x(1)^2 * x(2)^2
                    x(1) * x(2)^3                    
                                    ];
        
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

                                                        
    % *********************************************************************
    %
    % D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- MINIMAL BASIS FOR CRITIC, WITH SOME
    % ACTIVATION FUNCTIONS ADDED -- N_1 = 7 TERMS, ESSENTIAL ACTIVATION
    % FUNCTION x_1^2 ABSENT
    %             
    
    case 'vrabie_lewis_2009_hard_critic_7terms_no_x_1^2'  
        
        
        phix = [    
%                     x(1)^2
                    x(2)^2
                    x(2)^4                  % Opt. functs. end here  
                    x(1) * x(2) 
                    x(1)^4
                    x(1)^3 * x(2)
                    x(1)^2 * x(2)^2
                    x(1) * x(2)^3                    
                                    ];
        
        dphix = [   
%                     2 * x(1)            0
                    0                   2 * x(2)
                    0                   4 * x(2)^3  % Opt. functs. end here  
                    x(2)                x(1)  
                    4 * x(1)^3          0
                    3 * x(1)^2 * x(2)   x(1)^3
                    2 * x(1) * x(2)^2   2 * x(1)^2 * x(2)
                    x(2)^3              3 * x(1) * x(2)^2                    
                                                            ];                                                           
                                                        
                
    % *********************************************************************
    %
    % ALL METHOD COMPARISON OF D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- BASIS FOR CRITIC (DEBUGGING)
    %        
    
    case 'comp_vrabie_lewis_2009_hard_critic_debug'  
        
        
        phix = [    
                    x(1)^2
%                     x(1) * x(2)
                    x(2)^2
%                     x(1)^4
%                     x(1)^3 * x(2)
                    x(1)^2 * x(2)^2
%                     x(1) * x(2)^3
                    x(2)^4              
                                        ];
        
        dphix = [   
                    2 * x(1)            0
%                     x(2)                x(1)
                    0                   2 * x(2)
%                     4 * x(1)^3          0
%                     3 * x(1)^2 * x(2)   x(1)^3
                    2 * x(1) * x(2)^2   2 * x(1)^2 * x(2)
%                     x(2)^3              3 * x(1) * x(2)^2
                    0                   4 * x(2)^3          
                                                            ];        
                
                
    % *********************************************************************
    %
    % ALL METHOD COMPARISON OF D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- BASIS FOR ACTOR THAT DOESN'T USE g(x)
    %       
    
    case 'comp_vrabie_lewis_2009_hard_actor_no_g'
        
        phix = sin(x(1)) * [    x(2)
                                x(2)^3      ];
                 
%         phix = [    x(2)
%                     x(2)^3
%                     x(1)
%                     x(1)^2
%                     x(2)^2
%                     x(1) * x(2)
%                     x(1)^3     ];                
        
        dphix = [];     % Unused 
        
    % *********************************************************************
    %
    % ALL METHOD COMPARISON OF D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- BASIS FOR HAMILTONIAN NN FOR VI
    % ALGORITHM \Theta(x)
    %       
    
    case 'comp_vrabie_lewis_2009_hard_vi_theta'

        phix = sin(x(1))^2 * [      x(2)^2
                                    x(2)^4
                                    x(2)^6     ];        
                
        dphix = [];     % Unused        
           

    % *********************************************************************
    %
    % ALL METHOD COMPARISON OF D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- BASIS FOR HAMILTONIAN NN FOR VI
    % ALGORITHM \Theta(x) -- TERMS FOR Q(x) INCLUDED
    %       
    
    case 'comp_vrabie_lewis_2009_hard_vi_theta_w_Q'

            
        phix = [    
                    sin(x(1))^2 * x(2)^2
                    sin(x(1))^2 * x(2)^4
                    sin(x(1))^2 * x(2)^6 
                    x(1)^2
                    x(2)^2
                    x(2)^4
                                            ];
                
%         phix = [    
%                     sin(x(1))^2 * x(2)^2
%                     sin(x(1))^2 * x(2)^4
%                     sin(x(1))^2 * x(2)^6
%                     x(1)^2 + x(2)^2 + 2 * x(2)^4
%                                                     ];
                
        dphix = [];     % Unused          
        
        
    % *********************************************************************
    %
    % ALL METHOD COMPARISON OF D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- BASIS FOR HAMILTONIAN NN FOR VI
    % ALGORITHM \Psi(x)
    %             
    
    case 'comp_vrabie_lewis_2009_hard_vi_psi'
        
        phix = [    sin(x(1)) * x(2)
                    sin(x(1)) * x(2)^3      ];
        
        dphix = [];     % Unused          
        
    % *********************************************************************
    %
    % ALL METHOD COMPARISON OF D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2) -- BASIS FOR HAMILTONIAN NN FOR VI
    % ALGORITHM (TESTING)
    %            
    
    case 'comp_vrabie_lewis_2009_hard_vi_psi_test'
        
        phix = [    sin(x(1))^2 * x(2)^2
                    sin(x(1))^2 * x(2)^4
                    sin(x(1))^2 * x(2)^6
                    sin(x(1)) * x(2)
                    sin(x(1)) * x(2)^3      ];
        
        dphix = [];     % Unused         
        
    % *********************************************************************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %   
    
    otherwise
        
        error('*** ERROR: BASIS TAG NOT RECOGNIZED ***');

end