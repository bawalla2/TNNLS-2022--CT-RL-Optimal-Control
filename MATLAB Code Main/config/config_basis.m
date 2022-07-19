function basis = config_basis(basis, settings)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INITIALIZE BASIS PARAMETERS
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
% basis         (Struct) Contains basis parameters, partially initialized.
%               The fields of this input argument are algorithm-specific.
%               See respective algorithm m-file for details.
% settings      (Struct) Contains settings which are necessary to fully
%               initialize the basis. Has the fields:
%   sys         (Struct) Contains system parameters. See config_sys.m for
%               description of fields.
%   alg         (String) Algorithm tag. See config_main.m for options.
%   Q           (n x n matrix) State penalty function Q(x).
%   R           (m x m matrix) Control penalty matrix R = R^T > 0.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% basis         (Struct) Contains basis parameters, fully initialized.
%               The fields of this output argument are algorithm-specific.
%               See respective algorithm m-file for details.
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


% *************************************************************************
%
% UNPACK SETTINGS
% 
% *************************************************************************

% System used
sys = settings.sys;

% Algorithm used
alg = settings.alg;

% Cost structure
Q = settings.Q;
R = settings.R;


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% BASIS INITIALIZATION
%
% Get number of basis functions, and perform any other initialization
% necessary.
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
%
% CREATE CELL ARRAY OF BASIS OBJECTS
%
% Depending on the algorithm, there may be more than one basis used. Pack
% all bases into a cell array so that they may be manipulated/initialized
% one at a time.
% 
% *************************************************************************

switch alg

    % ***********************
    %
    % DEBUGGING
    %
    
    case 'debug'
        
        % Number of inputs to bases
        basis.numin = sys.n;

        basis_cell = {basis};
    
    % ***********************
    %
    % IRL
    %
    
    case 'irl'
        
        % Number of inputs to bases
        basis.numin = sys.n;
       
        basis_cell = {basis};

    % ***********************
    %
    % SPI
    %
    
    case 'spi'
  
        % Number of inputs to bases
        basis.numin = sys.n;
        
        basis_cell = {basis};
        
    % ***********************
    %
    % RADP (MATCHED UNCERTAINTY)
    %
    
    case 'radp_matched'
        
        % Number of inputs to bases
        basis.Phi.numin = sys.n;
        basis.Psi.numin = sys.n;
  
        basis_cell = {basis.Phi ; basis.Psi};       
        
        
    % ***********************
    %
    % RADP (UNMATCHED UNCERTAINTY)
    %
    
    case 'radp_unmatched'
        
        % Number of inputs to bases
        basis.Phi.numin = sys.n;
        basis.Psi.numin = sys.n;
        basis.Phi3.numin = sys.n + 1;   % Functions of (x,z)
        basis.Phi4.numin = sys.n;
  
        basis_cell = {basis.Phi ; basis.Psi ; basis.Phi3 ; basis.Phi4};


    % ***********************
    %
    % VI
    %
    
    case 'vi'        
        
        % Number of inputs to bases
        basis.Phi.numin = sys.n;
        basis.Psi.numin = sys.n;
        basis.Theta.numin = sys.n;
        
        basis_cell = {basis.Phi ; basis.Psi ; basis.Theta};

    % ***********************
    %
    % VI (TESTING)
    %
    
    case 'vi_test'
        
        % Number of inputs to bases
        basis.Phi.numin = sys.n;
        basis.Psi0.numin = sys.n;
        basis.Psi1.numin = sys.n;
        
        basis_cell = {basis.Phi ; basis.Psi0 ; basis.Psi1};
        
    % *********************************************************************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %   
    
    otherwise
        
        error('**** ERROR: ALGORITHM TAG NOT RECOGNIZED ***');  
       
end


% *************************************************************************
%
% INITIALIZE BASIS OBJECTS
% 
% *************************************************************************

% Number of bases used in this algorithm
numbasis = size(basis_cell, 1);

% Cell array to hold initialized basis objects
basis_cell_init = cell(numbasis, 1);

for i = 1:numbasis
   
    % Extract current basis
    basis_i = basis_cell{i};
    
    % *********************************************************************
    %
    % BASIS-SPECIFIC INITIALIZATION
    %
    % If the basis is of tag 'monomial' and type 'all', 'even', 'odd',
    % 'total_deg_K', the field 'dmat' needs to be initialized (see
    % description of eval_phi.m for details).
    %
    % *********************************************************************
    
%     % DEBUGGING
%     basis_i.tag = 'monomial';
%     basis_i.type = 'even';
%     basis_i.K = 4;
    
    if strcmp(basis_i.tag, 'monomial') 
        
        % Monomial basis used. See if 'dmat' needs to be initialized.
        
        if  strcmp(basis_i.type, 'all') ||...
            strcmp(basis_i.type, 'even') ||...
            strcmp(basis_i.type, 'odd') ||...
            strcmp(basis_i.type, 'total_deg_K')
        
            % 'dmat' needs to be initialized.
               
            % Max total degree of basis monomials
            K = basis_i.K;
            
            % Powers of values {0,1,...,K} to choose from for n factors of
            % the monomial
            a = repelem(0:K, sys.n);
            
            % All way to choose n elements from the vector a (before taking
            % permutations)
            d = nchoosek(a, sys.n);
            
            % Erase duplicates
            d = unique(d,'rows');
            
            % Get sum of columns along each row
            dsum = sum(d,2);
            
            % Remove rows which sum to zero and sum to > K
            d = d(((dsum > 0) & (dsum <= K)), :);
            
            % Get row sum of reduced matrix
            dsum = sum(d,2);
            
            % If only even/odd functions desired, remove odd/even
            % functions, respectively
            switch basis_i.type
                
                case 'even'
                    
                    % Remove odd-sum rows
                    d = d(mod(dsum,2) == 0, :);
                    
                case 'odd'
                    
                    % Remove even-sum rows
                    d = d(mod(dsum,2) == 1, :);
                    
                case 'total_deg_K'
                    
                    % Remove all rows of total degree ~= K
                    d = d(dsum == K, :);
                    
                case 'all'
                    
                    % Don't remove anything
                    
                % *********************************************************
                %
                % THROW ERROR IF TAG DOES NOT COME UP A MATCH
                %   
    
                otherwise
        
                    error(...
                     '*** ERROR: MONOMIAL BASIS TYPE NOT RECOGNIZED ***');  
                    
            end
            
            
            % Finally, take all permutations of the powers
            dmat = [];
            for j = 1:size(d,1)

                dmat = [    dmat
                            unique(perms(d(j,:)),'rows', 'stable')    ];

            end
            
            % Store dmat
            basis_i.dmat = dmat;
            
        else
            
            % The monomial basis is of type 'custom'. Thus, 'dmat' was
            % initialized manually already. Do nothing.
            
        end
        
        % *****************************************************************
        %
        % MONOMIAL BASIS WEIGHT INITIALIZATION
        %
        % Basis initialization options (see eval_phi.m for descriptions):
        %
        %   'none'  
        %   'custom'
        %   'lqr'  
        %     
        
        switch basis_i.IC_type
           
            case 'none'
                
                % Do nothing
                
            case 'custom'
                
                % Initialize the weights according to user settings
                dmat_ic = basis_i.dmat_ic;
                c0_ic = basis_i.c0_ic;
                c0 = config_basis_ICs(dmat, dmat_ic, c0_ic);
                
                % Store initial weights
                basis_i.c0 = c0;
                
            case 'lqr'
                
                % Initialize weights according to procedure outlined in
                % config_basis_ICs_lqr.m
                c0 = config_basis_ICs_lqr(sys, dmat, Q, R);
                
                % Store initial weights
                basis_i.c0 = c0;
                
            otherwise
                
                % *********************************************************
                %
                % THROW ERROR IF TAG DOES NOT COME UP A MATCH
                %
            
                error(['*** ERROR: MONOMIAL BASIS -- BASIS WEIGHT '...
                        'INITIALIZATION TAG NOT RECOGNIZED ***']);
                
        end
        
    end
    
    
    % *********************************************************************
    %
    % GET TOTAL NUMBER OF BASIS FUNCTIONS
    %
    % *********************************************************************
        
    
    % Evaluate activation functions at zero
    [phi0, ~] = eval_phi(zeros(basis_i.numin,1), basis_i);
    
    % Get number of activation functions
    basis_i.N = size(phi0, 1);
    

    % ***********************
    %
    % STORE INITIALIZED BASIS FUNCTION IN CELL ARRAY
    %
    
    basis_cell_init{i} = basis_i;
    
    
end


% *************************************************************************
%
% STORE INITIALIZED BASIS OBJECTS
%
% Now that each of the basis elements have been fully initialized, they
% need to be stored in the 'basis' struct. In effect, this we will take the
% basis objects stored in 'basis', which were partially initialized, and
% overwrite them with the corresponding objects in 'basis_cell_init', which
% have been fully initialized.
% 
% *************************************************************************


switch alg

    % ***********************
    %
    % DEBUGGING
    %
    
    case 'debug'   
        
        basis = basis_cell_init{1};
    
    % ***********************
    %
    % IRL
    %
    
    case 'irl'
       
        basis = basis_cell_init{1};

    % ***********************
    %
    % SPI
    %
    
    case 'spi'
  
        basis = basis_cell_init{1};
        
    % ***********************
    %
    % RADP (MATCHED UNCERTAINTY)
    %
    
    case 'radp_matched'
  
        basis.Phi = basis_cell_init{1};
        basis.Psi = basis_cell_init{2};  
            
    % ***********************
    %
    % RADP (UNMATCHED UNCERTAINTY)
    %
    
    case 'radp_unmatched'
  
        basis.Phi = basis_cell_init{1};
        basis.Psi = basis_cell_init{2}; 
        basis.Phi3 = basis_cell_init{3}; 
        basis.Phi4 = basis_cell_init{4}; 

    % ***********************
    %
    % VI
    %
    
    case 'vi'        
        
        basis.Phi = basis_cell_init{1};
        basis.Psi = basis_cell_init{2};
        basis.Theta = basis_cell_init{3};

    % ***********************
    %
    % VI (TESTING)
    %
    
    case 'vi_test'        
        
        basis.Phi = basis_cell_init{1};
        basis.Psi0 = basis_cell_init{2};
        basis.Psi1 = basis_cell_init{3};
        
    % *********************************************************************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %   
    
    otherwise
        
        error('**** ERROR: ALGORITHM TAG NOT RECOGNIZED ***');  
       
end


