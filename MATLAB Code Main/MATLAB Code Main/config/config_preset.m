function alg_settings = config_preset(preset, group_settings)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% SELECT ALGORITHM, SYSTEM, DESIGN PARAMETERS BASED ON PRESET
%
% Brent Wallace  
%
% 2021-11-06
%
% This program, given a desired example preset, handles all algorithm
% initialization/configuration.
%
% Algorithm options:
%
%   irl
%   spi
%   radp_matched
%   vi
%
% System options:
%
%   vrabie_lewis_2009_easy
%   vrabie_lewis_2009_hard          -- used in paper
%   vamvoudakis_lewis_2010_F16_lin
%   vamvoudakis_lewis_2010_nonlin
%   jiang_jiang_2014_engine
%   bian_jiang_2021_nonlin
%   cip                             -- used in paper  
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% alg_settings = config_preset(preset, group_settings)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% preset            (String) Algorithm/system preset for desired example.
% group_settings    (Struct) Contains system/design parameters to
%                   be shared across all presets in the desired group.
%                   E.g., if for this preset group all designs share the
%                   same Q, R matrices, those fields may be included in
%                   this struct. If not used, pass an empty array to this
%                   argument.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% alg_settings  (Struct) Algorithm settings/parameters for subsequent
%               execution according to desired preset (see respective
%               algorithm .m-file for details). 
%               
%               Regardless of the preset, alg_settings will have the
%               following fields:
%
%   plot_settings       (Struct) Contains plot settings for this preset.
%                       Has the following fields:
%       legend_entry    (String) A label for this specific preset to
%                       distinguish it from the other designs executed in
%                       the preset group. E.g., if a group executes
%                       different algorithms, legend_entry for one preset
%                       might be 'IRL'.
%       plotfolder      (String) Name of the folder to save plots to for
%                       this preset. This could be the preset tag, or any
%                       other convenient identifier.
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

% System
sys = group_settings.sys;


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CONFIGURE PRESET
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

switch preset

    
    %%    
    % *********************************************************************
    % *********************************************************************
    %
    % CASE STUDY I: 2ND ORDER ACADEMIC EXAMPLE, MINIMAL BASES
    %      
    % *********************************************************************
    % *********************************************************************
        
    % *********************************************************************
    %
    % CASE STUDY I: 2ND ORDER ACADEMIC EXAMPLE, MINIMAL BASES -- IRL
    %
    
    case 'CS_2ndorder_N1_3_irl'
    
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'irl';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
        
        % Number of policy iterations to perform in the PI algorithm before
        % terminating
        istar = group_settings.istar;
        
        % Number of simulations to run per iteration of the PI algorithm
        num_sims_per_iter = 1;      
        
        % Number of samples to collect per simulation for least squares
        % minimization at end of each iteration of PI algorithm
        l = 10;
        
        % Integral reinforcement interval length (sec)
%         tf = group_settings.tf;
        tf = 5;
        T = tf / (istar * num_sims_per_iter * ...
                                    l);        

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;
        
        % Probing noise
        noise.tag = '0';
%         noise = group_settings.noise;

        % Matrix of ICs. A total of 'istar' * 'num_sims_per_iter'
        % simulations are run in the algorithm. As many of the ICs as
        % desired can be set in this matrix. Else, ICs will be generated
        % manually.
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};        
        x0mat = x0';
        
        % Manner in which ICs are generated at the beginning of each new
        % simulation, after the algorithm has run out of ICs specified in
        % x0mat. Has the options:
        %       'rand'      Randomly generate new ICs. 
        %       'cont'      Continue ICs of new simulation as final values
        %                   of previous simulation.
        x0_behavior = 'cont';

        % ICs for critic NN
        c_0 = group_settings.basis_critic.c0;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
        
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'IRL';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY I: 2ND ORDER ACADEMIC EXAMPLE, MINIMAL BASES -- SPI
    %
    
    case 'CS_2ndorder_N1_3_spi'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'spi';   

        % Probing noise 
        noise = group_settings.noise;
        
        % Length of learning window [0, t_f].
        % I.e., time to insert probing noise for (sec).
%         tf = group_settings.tf;
        tf = 500;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;          
        
        % NN tuning gains
        alpha1 = 10;
        alpha2 = 10;
        
        % NN tuning parameters
        F1 = zeros(basis.N,1);   % Not used  
        F2 = 5*eye(basis.N);
        
        % System ICS
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end      
        x0 = group_settings.x0_cell{ind};           
        
        % ICs for critic NN
%         c_0 = zeros(basis.N,1);
%         c_0 = [1; 1; 0];
        c_0 = group_settings.basis_critic.c0;
        
        % ICs for actor NN
        w_0 = group_settings.basis_critic.c0;
                   
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'SPI';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

        
    % *********************************************************************
    %
    % CASE STUDY I: 2ND ORDER ACADEMIC EXAMPLE, MINIMAL BASES -- RADP
    % (MATCHED UNCERTAINTY)
    %        
        
    case 'CS_2ndorder_N1_3_radp_matched'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'radp_matched';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
                
        % Probing noise 
        noise = group_settings.noise;
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;         

        % Critic NN activation function parameters (cf. eqn. (13))
        Phi = group_settings.basis_critic; 
        
        % Actor NN activation function parameters (cf. eqn. (14))
        Psi = group_settings.basis_actor_no_g;
        
        % Pack basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;           
        
        % Setting to tell the algorithm if dynamic uncertainty is present
        % in the system.
        do_uncertainty = 0;
                
        % Number of samples to collect for least squares minimizations
        l = 50;
%         l = 10;
%         l = 5;
%         l = 100;
%         l = 1000;
        
        % Termination mode (see alg_radp_matched for details)
        term_mode = 'manual';

        % Tolerance \epsilon_1 to determine phase-1 termination (cf. eqn.
        % (37)) (if termination mode is convergence-based 'conv')
        eps1 = 1e-6;
        
        % Number of iterations i^* to terminate algorithm after (if
        % termination mode is 'manual')
        istar = group_settings.istar;
        
        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 15;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
%         tsim = 10;
             
        
        % Tag corresponding to the initial stabilizing policy u_0(x) (cf.
        % Assumption 4.2).
        u_0 = group_settings.u_0;
        
        % Robust redesign function \rho(s) (cf. under eqn. (31)).
        % NOTE: Not used (dynamic uncertainty not present).
        rho = '';
        
        % Initial conditions
        % IC for uncertain dynamics w (not used -- leave empty)
        w0 = [];         
        
        % IC for known dynamics x
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end     
        x0 = group_settings.x0_cell{ind};               
        
        % IC for actor NN w_0 \in R^{N_2}
%         w_0 = [-1 ; 1];  
        w_0 = 0 * ones(basis.Psi.N,1);
 
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'RADP';        

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY I: 2ND ORDER ACADEMIC EXAMPLE, MINIMAL BASES -- VI
    %
    
    case 'CS_2ndorder_N1_3_vi'
         
        % ***********************
        %
        % ALGORITHM AND SYSTEM
        %        
        
        alg = 'vi';
        
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %     
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R; 
        
        % Critic NN basis 
        Phi = group_settings.basis_critic;           
 
        % Hamiltonian NN basis \Psi(x)
        Psi = group_settings.basis_actor_no_g;        
        
        % Hamiltonian NN basis \Theta(x)
        Theta = group_settings.basis_hamiltonian.Theta;
                
        % Store basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;
        basis.Theta = Theta;
        
        % Probing noise tag
        noise = group_settings.noise;
        
        % Initial policy u_0(x)
%         u_0 = '0';                  % Not used
        u_0 = group_settings.u_0;
 
        % Length of time to tune weights for [0, s_f] (sec)
        sf = 125;
%         sf = 1000;

        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 10;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;        
                
        % Max stepsize for ode45
        maxstep = 1e-1;
%         maxstep = 1e-2;
%         maxstep = 1e-3;

        % Last function in Hamiltonian basis is r(x,u) (=1) or u^T R u (=0)
%         r1_R0 = 0;
        r1_R0 = 1;
        
        % Whether to include the last Hamiltonian activation function
        % explicitly in the basis
        include_in_Sigma = 0;
        
        % Whether to integrate Hamiltonian function (13) via ode45 (=1) or
        % manually (=0)
        int_H_ode45_1_man_0 = 0;
%         int_H_ode45_1_man_0 = 1;
        
        % Initial conditions x(0)
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end 
        x0 = group_settings.x0_cell{ind};
        
        % Initial critic weights c_0 \in R^{N_1}
%         c_0 = 0 * ones(Phi.N1,1);
        c_0 = [1 ;  1; 0];
%         c_0 = [0.5 ;  1; 1];

        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'VI';    
        
        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;        

    %%    
    % *********************************************************************
    % *********************************************************************
    %
    % CASE STUDY II: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 4
    % TERMS
    %      
    % *********************************************************************
    % *********************************************************************
        
    % *********************************************************************
    %
    % CASE STUDY II: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 4
    % TERMS -- IRL
    %
    
    case 'CS_2ndorder_N1_4_irl'
    
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'irl';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
        
        % Number of policy iterations to perform in the PI algorithm before
        % terminating
        istar = group_settings.istar;
        
        % Number of simulations to run per iteration of the PI algorithm
        num_sims_per_iter = 1;      
        
        % Number of samples to collect per simulation for least squares
        % minimization at end of each iteration of PI algorithm
        l = 10;
        
        % Integral reinforcement interval length (sec)
%         tf = group_settings.tf;
        tf = 5;
        T = tf / (istar * num_sims_per_iter * ...
                                    l);        

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;
        
        % Probing noise
        noise.tag = '0';
%         noise = group_settings.noise;

        % Matrix of ICs. A total of 'istar' * 'num_sims_per_iter'
        % simulations are run in the algorithm. As many of the ICs as
        % desired can be set in this matrix. Else, ICs will be generated
        % manually.
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};        
        x0mat = x0';
        
        % Manner in which ICs are generated at the beginning of each new
        % simulation, after the algorithm has run out of ICs specified in
        % x0mat. Has the options:
        %       'rand'      Randomly generate new ICs. 
        %       'cont'      Continue ICs of new simulation as final values
        %                   of previous simulation.
        x0_behavior = 'cont';

        % ICs for critic NN
        c_0 = group_settings.basis_critic.c0;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
        
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'IRL';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY II: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 4
    % TERMS -- SPI
    %
    
    case 'CS_2ndorder_N1_4_spi'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'spi';   

        % Probing noise 
        noise = group_settings.noise;
        
        % Length of learning window [0, t_f].
        % I.e., time to insert probing noise for (sec).
%         tf = group_settings.tf;
        tf = 500;
%         tf = 1;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;  
        
        % NN tuning gains
        alpha1 = 10;
        alpha2 = 10;
        
        % NN tuning parameters
        F1 = zeros(basis.N,1);   % Not used  
        F2 = 5*eye(basis.N);
        
        % System ICS
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end       
        x0 = group_settings.x0_cell{ind};           
        
        % ICs for critic NN
%         c_0 = zeros(basis.N,1);
%         c_0 = [1; 1; 0];
        c_0 = group_settings.basis_critic.c0;
        
        % ICs for actor NN
        w_0 = group_settings.basis_critic.c0;
                   
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'SPI';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

        
    % *********************************************************************
    %
    % CASE STUDY II: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 4
    % TERMS -- RADP (MATCHED UNCERTAINTY)
    %        
        
    case 'CS_2ndorder_N1_4_radp_matched'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'radp_matched';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
                
        % Probing noise 
        noise = group_settings.noise;
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;         

        % Critic NN activation function parameters (cf. eqn. (13))
        Phi = group_settings.basis_critic; 
        
        % Actor NN activation function parameters (cf. eqn. (14))
        Psi = group_settings.basis_actor_no_g;
        
        % Pack basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;           
        
        % Setting to tell the algorithm if dynamic uncertainty is present
        % in the system.
        do_uncertainty = 0;
                
        % Number of samples to collect for least squares minimizations
        l = 50;
%         l = 10;
%         l = 5;
%         l = 100;
%         l = 1000;
        
        % Termination mode (see alg_radp_matched for details)
        term_mode = 'manual';

        % Tolerance \epsilon_1 to determine phase-1 termination (cf. eqn.
        % (37)) (if termination mode is convergence-based 'conv')
        eps1 = 1e-6;
        
        % Number of iterations i^* to terminate algorithm after (if
        % termination mode is 'manual')
        istar = group_settings.istar;
        
        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 15;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
%         tsim = 10;
             
        
        % Tag corresponding to the initial stabilizing policy u_0(x) (cf.
        % Assumption 4.2).
        u_0 = group_settings.u_0;
        
        % Robust redesign function \rho(s) (cf. under eqn. (31)).
        % NOTE: Not used (dynamic uncertainty not present).
        rho = '';
        
        % Initial conditions
        % IC for uncertain dynamics w (not used -- leave empty)
        w0 = [];         
        
        % IC for known dynamics x
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end     
        x0 = group_settings.x0_cell{ind};               
        
        % IC for actor NN w_0 \in R^{N_2}
%         w_0 = [-1 ; 1];  
        w_0 = 0 * ones(basis.Psi.N,1);
 
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'RADP';        

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY II: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 4
    % TERMS -- VI
    %
    
    case 'CS_2ndorder_N1_4_vi'
         
        % ***********************
        %
        % ALGORITHM AND SYSTEM
        %        
        
        alg = 'vi';
        
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %     
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R; 
        
        % Critic NN basis 
        Phi = group_settings.basis_critic;           
 
        % Hamiltonian NN basis \Psi(x)
        Psi = group_settings.basis_actor_no_g;        
        
        % Hamiltonian NN basis \Theta(x)
        Theta = group_settings.basis_hamiltonian.Theta;
                
        % Store basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;
        basis.Theta = Theta;
        
        % Probing noise tag
        noise = group_settings.noise;
        
        % Initial policy u_0(x)
%         u_0 = '0';                  % Not used
        u_0 = group_settings.u_0;
 
        % Length of time to tune weights for [0, s_f] (sec)
        sf = 125;
%         sf = 1;

        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 5;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;        
        
        
        % Max stepsize for ode45
        maxstep = 1e-1;
%         maxstep = 1e-2;
%         maxstep = 1e-3;

        % Last function in Hamiltonian basis is r(x,u) (=1) or u^T R u (=0)
        r1_R0 = 1;
        
        % Whether to include the last Hamiltonian activation function
        % explicitly in the basis
        include_in_Sigma = 0;
 
        % Whether to integrate Hamiltonian function (13) via ode45 (=1) or
        % manually (=0)
        int_H_ode45_1_man_0 = 0;
%         int_H_ode45_1_man_0 = 1;        

        % Initial conditions x(0)
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};
        
        % Initial critic weights c_0 \in R^{N_1}
%         c_0 = 0 * ones(Phi.N,1);
        c_0 = [1 ;  1; zeros(basis.Phi.N - 2, 1)];
%         c_0 = [1 ;  1; 0; 0];
%         c_0 = [0.5 ;  1; 0.9; 0];
%         c_0 = [1 ;  1; 0];
%         c_0 = [0.5 ;  1; 1];

        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'VI';    
        
        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;                

        


    %%    
    % *********************************************************************
    % *********************************************************************
    %
    % CASE STUDY III: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_2^4 REMOVED
    %      
    % *********************************************************************
    % *********************************************************************
        
    % *********************************************************************
    %
    % CASE STUDY III: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_2^4 REMOVED -- IRL
    %
    
    case 'CS_2ndorder_N1_7_no_x2p4_irl'
    
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'irl';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
        
        % Number of policy iterations to perform in the PI algorithm before
        % terminating
        istar = group_settings.istar;
        
        % Number of simulations to run per iteration of the PI algorithm
        num_sims_per_iter = 1;      
        
        % Number of samples to collect per simulation for least squares
        % minimization at end of each iteration of PI algorithm
        l = 10;
%         l = 20;
        
        % Integral reinforcement interval length (sec)
%         tf = group_settings.tf;
        tf = 5;
        T = tf / (istar * num_sims_per_iter * ...
                                    l);      

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;
        
        % Probing noise
        noise.tag = '0';
%         noise = group_settings.noise;

        % Matrix of ICs. A total of 'istar' * 'num_sims_per_iter'
        % simulations are run in the algorithm. As many of the ICs as
        % desired can be set in this matrix. Else, ICs will be generated
        % manually.
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};        
        x0mat = x0';
        
        % Manner in which ICs are generated at the beginning of each new
        % simulation, after the algorithm has run out of ICs specified in
        % x0mat. Has the options:
        %       'rand'      Randomly generate new ICs. 
        %       'cont'      Continue ICs of new simulation as final values
        %                   of previous simulation.
        x0_behavior = 'cont';

        % ICs for critic NN
        c_0 = group_settings.basis_critic.c0;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim + (group_settings.tf - tf);
%         tsim = group_settings.tsim;
        
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'IRL';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY III: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_2^4 REMOVED -- SPI
    %
    
    case 'CS_2ndorder_N1_7_no_x2p4_spi'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'spi';   

        % Probing noise 
        noise = group_settings.noise;
        
        % Length of learning window [0, t_f].
        % I.e., time to insert probing noise for (sec).
%         tf = group_settings.tf;
%         tf = 37.75;
%         tf = 37.8;
%         tf = 25;
        tf = 100;
%         tf = 0.00001;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;  
        
        
        % NN tuning gains
        alpha1 = 10;
        alpha2 = 10;
        
        % NN tuning parameters
        F1 = zeros(basis.N,1);   % Not used  
        F2 = 5*eye(basis.N);
        
        % System ICS
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end       
        x0 = group_settings.x0_cell{ind};           
        
        % ICs for critic NN
%         c_0 = zeros(basis.N,1);
%         c_0 = [1; 1; 0];
        c_0 = group_settings.basis_critic.c0;
        
        % ICs for actor NN
        w_0 = group_settings.basis_critic.c0;
                   
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'SPI';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

        
    % *********************************************************************
    %
    % CASE STUDY III: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_2^4 REMOVED -- RADP (MATCHED
    % UNCERTAINTY)
    %        
        
    case 'CS_2ndorder_N1_7_no_x2p4_radp_matched'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'radp_matched';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
                
        % Probing noise 
        noise = group_settings.noise;
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;         

        % Critic NN activation function parameters (cf. eqn. (13))
        Phi = group_settings.basis_critic; 
        
        % Actor NN activation function parameters (cf. eqn. (14))
        Psi = group_settings.basis_actor_no_g;
        
        % Pack basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;   
        
        % Setting to tell the algorithm if dynamic uncertainty is present
        % in the system.
        do_uncertainty = 0;
                
        % Number of samples to collect for least squares minimizations
        l = 50;
%         l = 10;
%         l = 5;
%         l = 100;
%         l = 1000;
        
        % Termination mode (see alg_radp_matched for details)
        term_mode = 'manual';

        % Tolerance \epsilon_1 to determine phase-1 termination (cf. eqn.
        % (37)) (if termination mode is convergence-based 'conv')
        eps1 = 1e-6;
        
        % Number of iterations i^* to terminate algorithm after (if
        % termination mode is 'manual')
        istar = group_settings.istar;
        
        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 15;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
%         tsim = 10;
             
        
        % Tag corresponding to the initial stabilizing policy u_0(x) (cf.
        % Assumption 4.2).
        u_0 = group_settings.u_0;
        
        % Robust redesign function \rho(s) (cf. under eqn. (31)).
        % NOTE: Not used (dynamic uncertainty not present).
        rho = '';
        
        % Initial conditions
        % IC for uncertain dynamics w (not used -- leave empty)
        w0 = [];         
        
        % IC for known dynamics x
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end     
        x0 = group_settings.x0_cell{ind};               
        
        % IC for actor NN w_0 \in R^{N_2}
%         w_0 = [-1 ; 1];  
        w_0 = 0 * ones(basis.Psi.N,1);
 
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'RADP';        

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY III: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_2^4 REMOVED --  VI
    %
    
    case 'CS_2ndorder_N1_7_no_x2p4_vi'
         
        % ***********************
        %
        % ALGORITHM AND SYSTEM
        %        
        
        alg = 'vi';
        
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %     
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R; 
        
        % Critic NN basis 
        Phi = group_settings.basis_critic;           
 
        % Hamiltonian NN basis \Psi(x)
        Psi = group_settings.basis_actor_no_g;        
        
        % Hamiltonian NN basis \Theta(x)
        Theta = group_settings.basis_hamiltonian.Theta;
                
        % Store basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;
        basis.Theta = Theta;
        
        % Probing noise tag
        noise = group_settings.noise;
        
        % Initial policy u_0(x)
%         u_0 = '0';                  % Not used
        u_0 = group_settings.u_0;
 
        % Length of time to tune weights for [0, s_f] (sec)
%         sf = 50;
        sf = 0.01;

        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 5;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;        
        
        
        % Max stepsize for ode45
        maxstep = 1e-1;
%         maxstep = 1e-2;
%         maxstep = 1e-3;

        % Last function in Hamiltonian basis is r(x,u) (=1) or u^T R u (=0)
        r1_R0 = 1;
        
        % Whether to include the last Hamiltonian activation function
        % explicitly in the basis
        include_in_Sigma = 0;

        % Whether to integrate Hamiltonian function (13) via ode45 (=1) or
        % manually (=0)
        int_H_ode45_1_man_0 = 0;
%         int_H_ode45_1_man_0 = 1;        
        
        % Initial conditions x(0)
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};
        
        % Initial critic weights c_0 \in R^{N_1}
%         c_0 = 0 * ones(Phi.N,1);
        c_0 = [1 ;  1; zeros(basis.Phi.N - 2, 1)];
%         c_0 = [1 ;  1; 0; 0];
%         c_0 = [0.5 ;  1; 0.9; 0];
%         c_0 = [1 ;  1; 0];
%         c_0 = [0.5 ;  1; 1];

        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'VI';    
        
        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;                       

        
    %%    
    % *********************************************************************
    % *********************************************************************
    %
    % CASE STUDY IV: CART INVERTED PENDULUM SYSTEM
    %      
    % *********************************************************************
    % *********************************************************************
        
    % *********************************************************************
    %
    % CASE STUDY IV: CART INVERTED PENDULUM SYSTEM -- IRL
    %
    
    case 'CS_cip_irl'
    
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'irl';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
        
        % Number of policy iterations to perform in the PI algorithm before
        % terminating
        istar = group_settings.istar;
        
        % Number of simulations to run per iteration of the PI algorithm
        num_sims_per_iter = 1;      
        
        % Number of samples to collect per simulation for least squares
        % minimization at end of each iteration of PI algorithm
        l = 15;
%         l = 20;
%         l = 50;
        
        % Integral reinforcement interval length (sec)
%         tf = group_settings.tf;
        tf = 5;
        T = tf / (istar * num_sims_per_iter * ...
                                    l);      

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;
        
        % Probing noise
        noise.tag = '0';
%         noise = group_settings.noise;

        % Matrix of ICs. A total of 'istar' * 'num_sims_per_iter'
        % simulations are run in the algorithm. As many of the ICs as
        % desired can be set in this matrix. Else, ICs will be generated
        % manually.
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};        
        x0mat = x0';
        
        % Manner in which ICs are generated at the beginning of each new
        % simulation, after the algorithm has run out of ICs specified in
        % x0mat. Has the options:
        %       'rand'      Randomly generate new ICs. 
        %       'cont'      Continue ICs of new simulation as final values
        %                   of previous simulation.
        x0_behavior = 'cont';

        
        % ICs for critic NN
        c_0 = group_settings.basis_critic.c0; 
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
%         tsim = group_settings.tsim + (group_settings.tf - tf);
        tsim = group_settings.tsim;
        
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'IRL';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY IV: CART INVERTED PENDULUM SYSTEM -- SPI
    %
    
    case 'CS_cip_spi'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'spi';   

        % Probing noise 
        noise = group_settings.noise;
        
        % Length of learning window [0, t_f].
        % I.e., time to insert probing noise for (sec).
%         tf = group_settings.tf;
%         tf = 37.75;
%         tf = 37.8;
%         tf = 25;
        tf = 100;
%         tf = 0.00001;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic; 
        
        % NN tuning gains
%         alpha1 = 10;
        alpha1 = 1000;
        alpha2 = 10;
        
        % NN tuning parameters
        F1 = zeros(basis.N,1);   % Not used  
        F2 = 5*eye(basis.N);
        
        % System ICS
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end       
        x0 = group_settings.x0_cell{ind};           
        
        % ICs for critic NN
%         c_0 = zeros(basis.N,1);
%         c_0 = [1; 1; 0];
        c_0 = group_settings.basis_critic.c0;
        
        % ICs for actor NN
        w_0 = group_settings.basis_critic.c0;
                   
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'SPI';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

        
    % *********************************************************************
    %
    % CASE STUDY IV: CART INVERTED PENDULUM SYSTEM -- RADP (MATCHED
    % UNCERTAINTY)
    %        
        
    case 'CS_cip_radp_matched'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'radp_matched';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
                
        % Probing noise 
        noise = group_settings.noise;
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;         

        % Critic NN activation function parameters (cf. eqn. (13))
        Phi = group_settings.basis_critic; 
        
        % Actor NN activation function parameters (cf. eqn. (14))
        Psi = group_settings.basis_actor_no_g;
        
        % Pack basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;   
        
        % Setting to tell the algorithm if dynamic uncertainty is present
        % in the system.
        do_uncertainty = 0;
                
        % Number of samples to collect for least squares minimizations
%         l = 5;
%         l = 10;
        l = 50;
%         l = 75;
%         l = 100;
%         l = 1000;
        
        % Termination mode (see alg_radp_matched for details)
        term_mode = 'manual';

        % Tolerance \epsilon_1 to determine phase-1 termination (cf. eqn.
        % (37)) (if termination mode is convergence-based 'conv')
        eps1 = 1e-6;
        
        % Number of iterations i^* to terminate algorithm after (if
        % termination mode is 'manual')
        istar = group_settings.istar;
        
        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 15;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
%         tsim = 10;
             
        
        % Tag corresponding to the initial stabilizing policy u_0(x) (cf.
        % Assumption 4.2).
        u_0 = group_settings.u_0;
        
        % Robust redesign function \rho(s) (cf. under eqn. (31)).
        % NOTE: Not used (dynamic uncertainty not present).
        rho = '';
        
        % Initial conditions
        % IC for uncertain dynamics w (not used -- leave empty)
        w0 = [];         
        
        % IC for known dynamics x
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end     
        x0 = group_settings.x0_cell{ind};               
        
        % IC for actor NN w_0 \in R^{N_2}
%         w_0 = [-1 ; 1];  
        w_0 = 0 * ones(basis.Psi.N,1);
 
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'RADP';        

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY IV: CART INVERTED PENDULUM SYSTEM --  VI
    %
    
    case 'CS_cip_vi'
         
        % ***********************
        %
        % ALGORITHM AND SYSTEM
        %        
        
        alg = 'vi';
        
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %     
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R; 
        
        % Critic NN basis 
        Phi = group_settings.basis_critic;           
 
        % Hamiltonian NN basis \Psi(x)
        Psi = group_settings.basis_actor_no_g;        
        
        % Hamiltonian NN basis \Theta(x)
        Theta = group_settings.basis_hamiltonian.Theta;
                
        % Store basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;
        basis.Theta = Theta;
        
        % Probing noise tag
        noise = group_settings.noise;
        
        % Initial policy u_0(x)
%         u_0 = '0';                  % Not used
        u_0 = group_settings.u_0;
 
        % Length of time to tune weights for [0, s_f] (sec)
        sf = 1;
%         sf = 0.01;

        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 5;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;        
        
        
        % Max stepsize for ode45
        maxstep = 1e-1;
%         maxstep = 1e-2;
%         maxstep = 1e-3;

        % Last function in Hamiltonian basis is r(x,u) (=1) or u^T R u (=0)
        r1_R0 = 1;
        
        % Whether to include the last Hamiltonian activation function
        % explicitly in the basis
        include_in_Sigma = 0;

        % Whether to integrate Hamiltonian function (13) via ode45 (=1) or
        % manually (=0)
        int_H_ode45_1_man_0 = 0;
%         int_H_ode45_1_man_0 = 1;        
        
        % Initial conditions x(0)
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};
        
        % Initial critic weights c_0 \in R^{N_1}
%         c_0 = 0 * ones(Phi.N,1);
        c_0 = [1 ;  1; zeros(basis.Phi.N - 2, 1)];
%         c_0 = [1 ;  1; 0; 0];
%         c_0 = [0.5 ;  1; 0.9; 0];
%         c_0 = [1 ;  1; 0];
%         c_0 = [0.5 ;  1; 1];

        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'VI';    
        
        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;                               
        
        
    %%    
    % *********************************************************************
    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 8
    % TERMS
    %      
    % *********************************************************************
    % *********************************************************************
        
    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 8
    % TERMS -- IRL
    %
    
    case 'CS_2ndorder_N1_8_irl'
    
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'irl';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
        
        % Number of policy iterations to perform in the PI algorithm before
        % terminating
        istar = group_settings.istar;
        
        % Number of simulations to run per iteration of the PI algorithm
        num_sims_per_iter = 1;      
        
        % Number of samples to collect per simulation for least squares
        % minimization at end of each iteration of PI algorithm
        l = 10;
        
        % Integral reinforcement interval length (sec)
        tf = group_settings.tf;
        T = tf / (istar * num_sims_per_iter * ...
                                    l);        

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;
        
        % Probing noise
        noise.tag = '0';
%         noise = group_settings.noise;

        % Matrix of ICs. A total of 'istar' * 'num_sims_per_iter'
        % simulations are run in the algorithm. As many of the ICs as
        % desired can be set in this matrix. Else, ICs will be generated
        % manually.
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};        
        x0mat = x0';
        
        % Manner in which ICs are generated at the beginning of each new
        % simulation, after the algorithm has run out of ICs specified in
        % x0mat. Has the options:
        %       'rand'      Randomly generate new ICs. 
        %       'cont'      Continue ICs of new simulation as final values
        %                   of previous simulation.
        x0_behavior = 'cont';

        % ICs for critic NN
        c_0 = group_settings.basis_critic.c0;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
        
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'IRL';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 8
    % TERMS -- SPI
    %
    
    case 'CS_2ndorder_N1_8_spi'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'spi';   

        % Probing noise 
        noise = group_settings.noise;
        
        % Length of learning window [0, t_f].
        % I.e., time to insert probing noise for (sec).
%         tf = group_settings.tf;
        tf = 1500;
%         tf = 1;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;  
        
        % NN tuning gains
        alpha1 = 10;
        alpha2 = 10;
        
        % NN tuning parameters
        F1 = zeros(basis.N,1);   % Not used  
        F2 = 5*eye(basis.N);
        
        % System ICS
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end       
        x0 = group_settings.x0_cell{ind};           
        
        % ICs for critic NN
%         c_0 = zeros(basis.N,1);
%         c_0 = [1; 1; 0];
        c_0 = group_settings.basis_critic.c0;
        
        % ICs for actor NN
        w_0 = group_settings.basis_critic.c0;
                   
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'SPI';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

        
    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 8
    % TERMS -- RADP (MATCHED UNCERTAINTY)
    %        
        
    case 'CS_2ndorder_N1_8_radp_matched'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'radp_matched';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
                
        % Probing noise 
        noise = group_settings.noise;
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;         

        % Critic NN activation function parameters (cf. eqn. (13))
        Phi = group_settings.basis_critic; 
        
        % Actor NN activation function parameters (cf. eqn. (14))
        Psi = group_settings.basis_actor_no_g;
        
        % Pack basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;   
        
        % Setting to tell the algorithm if dynamic uncertainty is present
        % in the system.
        do_uncertainty = 0;
                
        % Number of samples to collect for least squares minimizations
        l = 50;
%         l = 10;
%         l = 5;
%         l = 100;
%         l = 1000;
        
        % Termination mode (see alg_radp_matched for details)
        term_mode = 'manual';

        % Tolerance \epsilon_1 to determine phase-1 termination (cf. eqn.
        % (37)) (if termination mode is convergence-based 'conv')
        eps1 = 1e-6;
        
        % Number of iterations i^* to terminate algorithm after (if
        % termination mode is 'manual')
        istar = group_settings.istar;
        
        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 15;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
%         tsim = 10;
             
        
        % Tag corresponding to the initial stabilizing policy u_0(x) (cf.
        % Assumption 4.2).
        u_0 = group_settings.u_0;
        
        % Robust redesign function \rho(s) (cf. under eqn. (31)).
        % NOTE: Not used (dynamic uncertainty not present).
        rho = '';
        
        % Initial conditions
        % IC for uncertain dynamics w (not used -- leave empty)
        w0 = [];         
        
        % IC for known dynamics x
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end     
        x0 = group_settings.x0_cell{ind};               
        
        % IC for actor NN w_0 \in R^{N_2}
%         w_0 = [-1 ; 1];  
        w_0 = 0 * ones(basis.Psi.N,1);
 
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'RADP';        

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 8
    % TERMS -- VI
    %
    
    case 'CS_2ndorder_N1_8_vi'
         
        % ***********************
        %
        % ALGORITHM AND SYSTEM
        %        
        
        alg = 'vi';
        
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %     
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R; 
        
        % Critic NN basis 
        Phi = group_settings.basis_critic;           
 
        % Hamiltonian NN basis \Psi(x)
        Psi = group_settings.basis_actor_no_g;        
        
        % Hamiltonian NN basis \Theta(x)
        Theta = group_settings.basis_hamiltonian.Theta;
                
        % Store basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;
        basis.Theta = Theta;
        
        % Probing noise tag
        noise = group_settings.noise;
        
        % Initial policy u_0(x)
%         u_0 = '0';                  % Not used
        u_0 = group_settings.u_0;
 
        % Length of time to tune weights for [0, s_f] (sec)
        sf = 50;
%         sf = 1;

        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 5;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;        
        
        
        % Max stepsize for ode45
        maxstep = 1e-1;
%         maxstep = 1e-2;
%         maxstep = 1e-3;

        % Last function in Hamiltonian basis is r(x,u) (=1) or u^T R u (=0)
        r1_R0 = 1;
        
        % Whether to include the last Hamiltonian activation function
        % explicitly in the basis
        include_in_Sigma = 0;

        % Whether to integrate Hamiltonian function (13) via ode45 (=1) or
        % manually (=0)
        int_H_ode45_1_man_0 = 0;
%         int_H_ode45_1_man_0 = 1;        
        

        % Initial conditions x(0)
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};
        
        % Initial critic weights c_0 \in R^{N_1}
%         c_0 = 0 * ones(Phi.N,1);
        c_0 = [1 ;  1; zeros(basis.Phi.N - 2, 1)];
%         c_0 = [1 ;  1; 0; 0];
%         c_0 = [0.5 ;  1; 0.9; 0];
%         c_0 = [1 ;  1; 0];
%         c_0 = [0.5 ;  1; 1];

        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'VI';    
        
        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;                       
                
        

    %%    
    % *********************************************************************
    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_1^2 REMOVED
    %      
    % *********************************************************************
    % *********************************************************************
        
    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_1^2 REMOVED -- IRL
    %
    
    case 'CS_2ndorder_N1_7_no_x1p2_irl'
    
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'irl';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
        
        % Number of policy iterations to perform in the PI algorithm before
        % terminating
        istar = group_settings.istar;
        
        % Number of simulations to run per iteration of the PI algorithm
        num_sims_per_iter = 1;      
        
        % Number of samples to collect per simulation for least squares
        % minimization at end of each iteration of PI algorithm
        l = 10;
        
        % Integral reinforcement interval length (sec)
        tf = group_settings.tf;
        T = tf / (istar * num_sims_per_iter * ...
                                    l);        

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;
        
        % Probing noise
        noise.tag = '0';
%         noise = group_settings.noise;

        % Matrix of ICs. A total of 'istar' * 'num_sims_per_iter'
        % simulations are run in the algorithm. As many of the ICs as
        % desired can be set in this matrix. Else, ICs will be generated
        % manually.
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};        
        x0mat = x0';
        
        % Manner in which ICs are generated at the beginning of each new
        % simulation, after the algorithm has run out of ICs specified in
        % x0mat. Has the options:
        %       'rand'      Randomly generate new ICs. 
        %       'cont'      Continue ICs of new simulation as final values
        %                   of previous simulation.
%         x0_behavior = 'rand';
        x0_behavior = 'cont';

        % ICs for critic NN
        c_0 = group_settings.basis_critic.c0;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
        
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'IRL';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_1^2 REMOVED -- SPI
    %
    
    case 'CS_2ndorder_N1_7_no_x1p2_spi'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'spi';   

        % Probing noise 
        noise = group_settings.noise;
        
        % Length of learning window [0, t_f].
        % I.e., time to insert probing noise for (sec).
%         tf = group_settings.tf;
        tf = 750;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;

        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;
        
        % Basis
        basis = group_settings.basis_critic;  
        
        % NN tuning gains
        alpha1 = 10;
        alpha2 = 10;
        
        % NN tuning parameters
        F1 = zeros(basis.N,1);   % Not used  
        F2 = 5*eye(basis.N);
        
        % System ICS
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end       
        x0 = group_settings.x0_cell{ind};           
        
        % ICs for critic NN
%         c_0 = zeros(basis.N,1);
%         c_0 = [1; 1; 0];
        c_0 = group_settings.basis_critic.c0;
        
        % ICs for actor NN
        w_0 = group_settings.basis_critic.c0;
                   
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'SPI';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

        
    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_1^2 REMOVED -- RADP (MATCHED
    % UNCERTAINTY)
    %        
        
    case 'CS_2ndorder_N1_7_no_x1p2_radp_matched'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'radp_matched';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
                
        % Probing noise 
        noise = group_settings.noise;
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;         

        % Critic NN activation function parameters (cf. eqn. (13))
        Phi = group_settings.basis_critic; 
        
        % Actor NN activation function parameters (cf. eqn. (14))
        Psi = group_settings.basis_actor_no_g;
        
        % Pack basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;   
        
        % Setting to tell the algorithm if dynamic uncertainty is present
        % in the system.
        do_uncertainty = 0;
                
        % Number of samples to collect for least squares minimizations
        l = 50;
%         l = 10;
%         l = 5;
%         l = 100;
%         l = 1000;
        
        % Termination mode (see alg_radp_matched for details)
        term_mode = 'manual';

        % Tolerance \epsilon_1 to determine phase-1 termination (cf. eqn.
        % (37)) (if termination mode is convergence-based 'conv')
        eps1 = 1e-6;
        
        % Number of iterations i^* to terminate algorithm after (if
        % termination mode is 'manual')
        istar = group_settings.istar;
        
        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 15;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
%         tsim = 10;
             
        
        % Tag corresponding to the initial stabilizing policy u_0(x) (cf.
        % Assumption 4.2).
        u_0 = group_settings.u_0;
        
        % Robust redesign function \rho(s) (cf. under eqn. (31)).
        % NOTE: Not used (dynamic uncertainty not present).
        rho = '';
        
        % Initial conditions
        % IC for uncertain dynamics w (not used -- leave empty)
        w0 = [];         
        
        % IC for known dynamics x
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end     
        x0 = group_settings.x0_cell{ind};               
        
        % IC for actor NN w_0 \in R^{N_2}
%         w_0 = [-1 ; 1];  
        w_0 = 0 * ones(basis.Psi.N,1);
 
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'RADP';        

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;           

    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_1^2 REMOVED --  VI
    %
    
    case 'CS_2ndorder_N1_7_no_x1p2_vi'
         
        % ***********************
        %
        % ALGORITHM AND SYSTEM
        %        
        
        alg = 'vi';
        
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %     
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R; 
        
        % Critic NN basis 
        Phi = group_settings.basis_critic;           
 
        % Hamiltonian NN basis \Psi(x)
        Psi = group_settings.basis_actor_no_g;        
        
        % Hamiltonian NN basis \Theta(x)
        Theta = group_settings.basis_hamiltonian.Theta;
                
        % Store basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;
        basis.Theta = Theta;
        
        % Probing noise tag
        noise = group_settings.noise;
        
        % Initial policy u_0(x)
%         u_0 = '0';                  % Not used
        u_0 = group_settings.u_0;
 
        % Length of time to tune weights for [0, s_f] (sec)
        sf = 50;

        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 5;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;        
        
        
        % Max stepsize for ode45
        maxstep = 1e-1;
%         maxstep = 1e-2;
%         maxstep = 1e-3;

        % Last function in Hamiltonian basis is r(x,u) (=1) or u^T R u (=0)
        r1_R0 = 1;
        
        % Whether to include the last Hamiltonian activation function
        % explicitly in the basis
        include_in_Sigma = 0;

        % Whether to integrate Hamiltonian function (13) via ode45 (=1) or
        % manually (=0)
        int_H_ode45_1_man_0 = 0;
%         int_H_ode45_1_man_0 = 1;        

        % Initial conditions x(0)
        ind = mod(group_settings.presetcount, group_settings.numICs); 
        if ind == 0
            ind = group_settings.numICs;
        end
        x0 = group_settings.x0_cell{ind};
        
        % Initial critic weights c_0 \in R^{N_1}
%         c_0 = 0 * ones(Phi.N,1);
        c_0 = [1 ;  1; zeros(basis.Phi.N - 2, 1)];
%         c_0 = [1 ;  1; 0; 0];
%         c_0 = [0.5 ;  1; 0.9; 0];
%         c_0 = [1 ;  1; 0];
%         c_0 = [0.5 ;  1; 1];

        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = 'VI';    
        
        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;                       
                        
        
        
    
    %%
    % *********************************************************************
    %
    % VRABIE, LEWIS (2009) -- NONLINEAR "HARD" EXAMPLE
    %
    
    case 'vrabie_lewis_2009_hard_ex'
        
        %%   
        % ***********************
        %
        % ALGORITHM 
        %        
        
        alg = 'irl';

        % State penalty function
        Q = 'vrabie_lewis_2009_hard_ex';
        
        % Control penalty matrix
        R = 1;          
        
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %    
   
        
        % Number of policy iterations to perform in the PI algorithm before
        % terminating
        istar = 5;
        
        % Number of simulations to run per iteration of the PI algorithm
        num_sims_per_iter = 8;      
        
        % Number of samples to collect per simulation for least squares
        % minimization at end of each iteration of PI algorithm
        l = 5;
        
        
%         % DEBUGGING: Compare to Vrabie, Lewis (2009) code settings
%         istar = 5;
%         num_sims_per_iter = 4;
%         l = 5;
     
%         istar = 5;
%         num_sims_per_iter = 1;
%         l = 10;

        % Integral reinforcement interval length (sec)
        tf = 20;
        T = tf / (istar * num_sims_per_iter * ...
                                    l); 

        % Matrix of ICs. A total of 'istar' * 'num_sims_per_iter'
        % simulations are run in the algorithm. As many of the ICs as
        % desired can be set in this matrix. Else, ICs will be generated
        % manually.
        x0mat = [];

        % Manner in which ICs are generated at the beginning of each new
        % simulation, after the algorithm has run out of ICs specified in
        % x0mat. Has the options:
        %       'rand'      Randomly generate new ICs. 
        %       'cont'      Continue ICs of new simulation as final values
        %                   of previous simulation.
        x0_behavior = 'rand';        
        
        % ICs for critic NN
        c_0 = [0 ; 0 ; 3/2 ; 0 ; 0 ; -0.1 ; 0 ; 3];
        
        % Basis
        basis.tag = 'order_2_degree_4';
        
        % Initialize the basis parameters
        b_sett.sys = sys;
        b_sett.alg = alg;
        b_sett.Q = Q;
        b_sett.R = R;
        basis = config_basis(basis, b_sett);      

        % Probing noise
        noise.tag = '0';        % No probing noise used

        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = 10;         
        
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = '';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = preset;        
        
    % *********************************************************************
    %
    % VAMVOUDAKIS, LEWIS (2010) -- LINEAR F16 EXAMPLE
    %
    
    case 'vamvoudakis_lewis_2010_F16_lin_ex'
        
        %%
        % ***********************
        %
        % ALGORITHM AND SYSTEM
        %        
        
        alg = 'spi';       

        % Probing noise tag
        noise.tag = 'vamvoudakis_lewis_2010_ex';
        
        % Length of learning window [0, t_f].
        % I.e., time to insert probing noise for (sec).
        tf = 750;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = 50;        

        % State penalty matrix
        Q = eye(3);
        
        % Control penalty matrix
        R = 1;
        
        % Basis
        basis.tag = 'order_3_degree_2';
        
        % Initialize the basis parameters
        b_sett.sys = sys;
        b_sett.alg = alg;
        b_sett.Q = Q;
        b_sett.R = R;        
        basis = config_basis(basis, b_sett);
        
        % NN tuning gains
        alpha1 = 50;
        alpha2 = 1;
        
        % NN tuning parameters
        F1 = -5*ones(basis.N,1);     
        F2 = 10*eye(basis.N);
        
        % Initial conditions
        c_0 = ones(basis.N,1);
        w_0 = rand(basis.N,1);
        x0 = [1; -1; 1];
                
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = '';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = preset;        
        
    % *********************************************************************
    %
    % VAMVOUDAKIS, LEWIS (2010) -- NONLINEAR EXAMPLE
    %
    
    case 'vamvoudakis_lewis_2010_nonlin_ex'  
   
        %%
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'spi';       

        % Probing noise tag
        noise.tag = 'vamvoudakis_lewis_2010_ex_nonlin';

        % Length of learning window [0, t_f].
        % I.e., time to insert probing noise for (sec).
        tf = 80;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = 20;        
        
        % State penalty matrix
        Q = eye(2);
        
        % Control penalty matrix
        R = 1;
        
        % Basis
        basis.tag = 'order_2_degree_2';
        
        % Initialize the basis parameters
        b_sett.sys = sys;
        b_sett.alg = alg;
        b_sett.Q = Q;
        b_sett.R = R;
        basis = config_basis(basis, b_sett);
        
        % NN tuning gains
        alpha1 = 5;
        alpha2 = 1;
        
        % NN tuning parameters
        F1 = 1*ones(basis.N,1);     
        F2 = 1*eye(basis.N);
                
        % Initial conditions
        c_0 = ones(basis.N,1);
        w_0 = ones(basis.N,1);
        x0 = [1; 1];

%         % DEBUGGING: Trying to get original Wa update law to stabilize
%         alpha1 = 10;
%         alpha2 = 5;
%         F1 = 10*ones(basis.dim,1);     
%         F2 = 100*eye(basis.dim);
%         c_0 = [0.5; 0; 1];
%         w_0 = [0.5; 0; 1];   
%         tf = 10;

        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = '';

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = preset;        
        
    % *********************************************************************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE EXAMPLE
    %        
        
    case 'jiang_jiang_2014_engine_ex'
        
        %%
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'radp_unmatched';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
                
        % Probing noise tag
        noise.tag = 'jiang_jiang_2014_engine_ex';
        
        % Number of samples to collect for least squares minimizations
%         l = 50;
%         l = 10;
%         l = 5;
        l = 100;
%         l = 1000;
        
        % Tolerance \epsilon_1 to determine phase-1 termination (cf. eqn.
        % (59))
        eps1 = 1e-6;

        % Online data collection window [0, t_f]
        tf = 10;   
%         tf = 1000;
        
        % Time to simulate for after getting final robust policy
        tsim = 30;
%         tsim = 10;
%         tsim = 0.001;

        % State penalty function tag
        Q = 'jiang_jiang_2014_engine_ex';
        
        % Constant \epsilon > 0 such that Q(x) - \epsilon^2 ||x||^2 is a
        % positive-definite function (cf. Assumption 3.3, Sec. III. B.)
        % NOTE: In this case Q(x) = 4(x^2 + x^3 + x^4), so \epsilon \in
        % (0,3] works
        eps = sqrt(0.1);
        
        % Control penalty gain
        R = 1;        

        % Critic NN activation function parameters (cf. eqn. (13))
        Phi.tag = 'jiang_jiang_2014_engine_ex_Phi1';
        
        % Actor NN activation function parameters (cf. eqn. (14))
        Psi.tag = 'jiang_jiang_2014_engine_ex_Phi2';
        
        % NN activation function parameters for approximation of
        % \hat{f}_1(x,z) (cf. eqn. (43))
        Phi3.tag = 'jiang_jiang_2014_engine_ex_Phi3';
        
        % NN activation function parameters for approximation of
        % \hat{g}_1(x) (cf. eqn. (44))
        Phi4.tag = 'jiang_jiang_2014_engine_ex_Phi4';
        
        % Pack basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;
        basis.Phi3 = Phi3;
        basis.Phi4 = Phi4;
        
        % Initialize the basis parameters
        b_sett.sys = sys;
        b_sett.alg = alg;
        b_sett.Q = Q;
        b_sett.R = R;
        basis = config_basis(basis, b_sett);
        
        % Tag corresponding to the initial stabilizing policy u_0(x) (cf.
        % Assumption 4.2).
        u_0 = 'jiang_jiang_2014_engine_ex';
        
        % Robust redesign function \rho(s) (cf. under eqn. (31)).
        rho = 'jiang_jiang_2014_engine_ex';
        
        % Initial conditions
        w0 = 1;                 % IC for uncertain dynamics w
        x0 = 2;                 % IC for known dynamics x
        z0 = -0.1;              % IC for known dynamics z
        
        w_0 = [1; zeros(basis.Psi.N - 1, 1)];       % IC for actor NN
 
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = '';        

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = preset;        
        
    % *********************************************************************
    %
    % T. BIAN, Z.-P. JIANG (2021) -- NONLINEAR EXAMPLE
    %
    
    case 'bian_jiang_2021_nonlin'
        
        %%
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'vi'; 
        
        % Order in which the control signal enters the system
        sys.upow = 3;
        
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %     
        
        % State, control penalty functions
        Q = 'bian_jiang_2021_nonlin_ex';
        R = 'bian_jiang_2021_nonlin_ex';
        
        % Critic NN basis 
        Phi.tag = 'order_2_degree_2';
        
        % Hamiltonian NN basis \Psi(x)
        Psi.tag = 'bian_jiang_2021_nonlin_psi';
        
        % Hamiltonian NN basis \Theta(x)
        Theta.tag = 'bian_jiang_2021_nonlin_theta';
           
        % Store basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;
        basis.Theta = Theta;
        
        % Initialize the basis parameters
        b_sett.sys = sys;
        b_sett.alg = alg;
        b_sett.Q = Q;
        b_sett.R = R;
        basis = config_basis(basis, b_sett);
        
        % Probing noise tag
        noise.tag = 'bian_jiang_2021_nonlin';
        
        % Initial policy u_0(x) (not used)
        u_0.tag = '0';
 
        % Weight learning time s_f (sec)
        sf = 100;
        
        % Length of learning window [0, t_f] (sec)
        tf = 1.8;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = 10 - tf;                  
        
        % Forward Euler approximation for NN weight updates (cf. Sec. IV.
        % D.)
        h_k = 0.003;
%         h_k = 1e-2;
        
        % Max stepsize for ode45
        maxstep = 1e-3;

        % Last function in Hamiltonian basis is r(x,u) (=1) or u^T R u (=0)
        r1_R0 = 0;

        % Whether to integrate Hamiltonian function (13) via ode45 (=1) or
        % manually (=0)
        int_H_ode45_1_man_0 = 0;
%         int_H_ode45_1_man_0 = 1;        
        
        % Whether to include the last Hamiltonian activation function
        % explicitly in the basis
        include_in_Sigma = 0;
        
        
        % Initial conditions x(0)
        x0 = [-2.9; -2.9];
        
        % Initial critic weights c_0 \in R^{N_1}
        c_0 = zeros(3,1);

        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry (not used)
        legend_entry = '';    
        
        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = preset;

    %%    
    % *********************************************************************
    % *********************************************************************
    %
    % RADP ON 2ND ORDER SYSTEM -- SWEEP x_0
    %      
    % *********************************************************************
    % *********************************************************************        
          
        
    case 'radp_2ndorder_sweep_x0'
        
        
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'radp_matched';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
                
        % Probing noise 
        noise = group_settings.noise;
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;         

        % Critic NN activation function parameters (cf. eqn. (13))
        Phi = group_settings.basis_critic; 
        
        % Actor NN activation function parameters (cf. eqn. (14))
        Psi = group_settings.basis_actor_no_g;
        
        % Pack basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;   
        
        % Setting to tell the algorithm if dynamic uncertainty is present
        % in the system.
        do_uncertainty = 0;
                
        % Number of samples to collect for least squares minimizations
%         l = 50;
%         l = 10;
%         l = 5;
        l = 100;
%         l = 1000;
        
        % Termination mode
        term_mode = 'manual';
        
        % Max number of iterations until convergence
        istar = group_settings.istar;

        % Tolerance \epsilon_1 to determine phase-1 termination (cf. eqn.
        % (37))
        eps1 = 1e-6;
        
        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
%         tf = 10;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
%         tsim = 0.1;
             
        
        % Tag corresponding to the initial stabilizing policy u_0(x) (cf.
        % Assumption 4.2).
        u_0 = group_settings.u_0;
        
        % Robust redesign function \rho(s) (cf. under eqn. (31)).
        % NOTE: Not used (dynamic uncertainty not present).
        rho = '';
        
        % Initial conditions
        % IC for uncertain dynamics w (not used -- leave empty)
        w0 = [];      
        
        % IC for known dynamics x
        presetcount = group_settings.presetcount;
        x0 = group_settings.x0_cell{presetcount};                
        
        % IC for actor NN w_0 \in R^{N_2}
%         w_0 = [-1 ; 1];  
        w_0 = 0 * ones(basis.Psi.N,1);
 
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'RADP';        

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;                 

        
    %%    
    % *********************************************************************
    % *********************************************************************
    %
    % RADP ON CART INVERTED PENDULUM SYSTEM -- SWEEP x_0
    %      
    % *********************************************************************
    % *********************************************************************        
          
        
    case 'radp_cip_sweep_x0'
                
        % ***********************
        %
        % ALGORITHM
        %        
        
        alg = 'radp_matched';
  
        % ***********************
        %
        % SETTINGS AND DESIGN PARAMETERS
        %       
                
        % Probing noise
        noise = group_settings.noise;
        
        % State penalty function
        Q = group_settings.Q;
        
        % Control penalty matrix
        R = group_settings.R;         

        % Critic NN activation function parameters (cf. eqn. (13))
        Phi = group_settings.basis_critic; 
        
        % Actor NN activation function parameters (cf. eqn. (14))
        Psi = group_settings.basis_actor_no_g;
        
        % Pack basis parameters
        basis.Phi = Phi;
        basis.Psi = Psi;   
        
        % Setting to tell the algorithm if dynamic uncertainty is present
        % in the system.
        do_uncertainty = 0;
                
        % Number of samples to collect for least squares minimizations
%         l = 50;
%         l = 10;
%         l = 5;
        l = 100;
%         l = 1000;
        
        % Tolerance \epsilon_1 to determine phase-1 termination (cf. eqn.
        % (37))
        eps1 = 1e-6;
        
        % Length of learning window [0, t_f] (sec)
        tf = group_settings.tf;
        
        % Time to simululate for after learning (i.e, [t_f, t_f + t_sim])
        % (sec)
        tsim = group_settings.tsim;
                  
        % Tag corresponding to the initial stabilizing policy u_0(x) (cf.
        % Assumption 4.2).
        u_0 = group_settings.u_0;
        
        % Robust redesign function \rho(s) (cf. under eqn. (31)).
        % NOTE: Not used (dynamic uncertainty not present).
        rho = '';
        
        % Initial conditions
        % IC for uncertain dynamics w (not used -- leave empty)
        w0 = [];      
        
        % IC for known dynamics x
        presetcount = group_settings.presetcount;
        xtheta0 = group_settings.x0_cell{presetcount};
        x0 = [xtheta0(1) ; 0 ; xtheta0(2) ; 0];                
        
        % IC for actor NN w_0 \in R^{N_2}
        w_0 = 0 * ones(basis.Psi.N,1);
 
        % ***********************
        %
        % PLOT SETTINGS
        %    
        
        % Legend entry
        legend_entry = 'RADP';        

        % Plot folder name. Can use preset tag, or anything else.
        plotfolder = legend_entry;                 
                
        
        


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





%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% STORE ALGORITHM SETTINGS/DESIGN PARAMETERS
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
%
% GENERAL SETTINGS
%
% *************************************************************************

% Preset tag
alg_settings.preset = preset;

% Plot settings -- general
plot_settings.legend_entry = legend_entry;
plot_settings.plotfolder = plotfolder;
% plot_settings.sys_settings = sys_settings;

% Write plot settings
alg_settings.plot_settings = plot_settings;


% *************************************************************************
%
% ALGORITHM-SPECIFIC SETTINGS
%
% *************************************************************************

switch alg

    % *********************************************************************
    %
    % IRL
    %
    
    case 'irl'
       
        alg_settings.sys = sys;
        alg_settings.alg = alg;
 
        alg_settings.Q = Q;
        alg_settings.R = R;
        alg_settings.basis = basis;
        alg_settings.tf = tf;
        alg_settings.noise = noise;
        
        alg_settings.T = T;
        alg_settings.istar = istar;
        alg_settings.num_sims_per_iter = num_sims_per_iter;
        alg_settings.l = l;
        
        alg_settings.x0mat = x0mat;
        alg_settings.x0_behavior = x0_behavior;
        alg_settings.c_0 = c_0;
        
        alg_settings.tsim = tsim;

    % *********************************************************************
    %
    % SPI
    %
    
    case 'spi'
  
        alg_settings.sys = sys;
        alg_settings.alg = alg;
              
        alg_settings.Q = Q;
        alg_settings.R = R;
        alg_settings.basis = basis;
        alg_settings.noise = noise;
        alg_settings.tf = tf;
        alg_settings.tsim = tsim;
        
        alg_settings.alpha1 = alpha1;
        alg_settings.alpha2 = alpha2;
        alg_settings.F1 = F1;
        alg_settings.F2 = F2;
        
        alg_settings.c_0 = c_0;
        alg_settings.w_0 = w_0;
        alg_settings.x0 = x0;
 
    % *********************************************************************
    %
    % RADP (MATCHED UNCERTAINTY)
    %
    
    case 'radp_matched'
  
        alg_settings.sys = sys;
        alg_settings.alg = alg;
        
        alg_settings.Q = Q;
        alg_settings.R = R;
        alg_settings.basis = basis;
        alg_settings.noise = noise;
        alg_settings.do_uncertainty = do_uncertainty;
        alg_settings.rho = rho;
        
        alg_settings.l = l;
        
        alg_settings.term_mode = term_mode;
        switch term_mode
            case 'conv'
                alg_settings.eps1 = eps1;
            case 'manual'
                alg_settings.istar = istar;
        end        
        
        
        alg_settings.tf = tf;
        alg_settings.tsim = tsim;
        
        alg_settings.u_0 = u_0;
        
        alg_settings.w0 = w0;
        alg_settings.x0 = x0;  
        
        alg_settings.w_0 = w_0;        
        
        
    % *********************************************************************
    %
    % RADP (UNMATCHED UNCERTAINTY)
    %
    
    case 'radp_unmatched'
  
        alg_settings.sys = sys;
        alg_settings.alg = alg;
        
        alg_settings.Q = Q;
        alg_settings.R = R;
        alg_settings.basis = basis;
        alg_settings.noise = noise;
        alg_settings.rho = rho;
        
        alg_settings.l = l;
        alg_settings.eps = eps;
        
        alg_settings.term_mode = term_mode;
        switch term_mode
            case 'conv'
                alg_settings.eps1 = eps1;
            case 'manual'
                alg_settings.istar = istar;
        end

        alg_settings.tf = tf;
        alg_settings.tsim = tsim;
        
        alg_settings.u_0 = u_0;
        
        alg_settings.w0 = w0;
        alg_settings.x0 = x0;
        alg_settings.z0 = z0;    
        
        alg_settings.w_0 = w_0;


    % *********************************************************************
    %
    % VI
    %
    
    case 'vi'        
        
        alg_settings.sys = sys;
        alg_settings.alg = alg;        

        alg_settings.Q = Q;
        alg_settings.R = R;
        alg_settings.basis = basis;
        alg_settings.noise = noise;
        alg_settings.u_0 = u_0;
        
        alg_settings.sf = sf;
        alg_settings.tf = tf;
        alg_settings.tsim = tsim;
        alg_settings.maxstep = maxstep;
        alg_settings.r1_R0 = r1_R0;
        alg_settings.include_in_Sigma = include_in_Sigma;
        alg_settings.int_H_ode45_1_man_0 = int_H_ode45_1_man_0;

        alg_settings.x0 = x0;
        alg_settings.c_0 = c_0;


        
    % *********************************************************************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %   
    
    otherwise
        
        error('**** ERROR: ALGORITHM TAG NOT RECOGNIZED ***');  
       
end
