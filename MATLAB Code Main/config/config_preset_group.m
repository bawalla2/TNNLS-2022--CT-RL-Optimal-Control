function [preset_list, group_settings] = config_preset_group(preset_group)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CONFIGURE PRESET GROUP
%
% Brent Wallace  
%
% 2022-02-16
%
% This program, given a desired preset group, initializes each of the
% presets to run in the group. See main.m for a list of preset groups.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% [preset_list, group_settings] = config_preset_group(preset_group)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% preset_group      (String) Tag corresponding to the preset group to run.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% preset_list       ('numpresets' x 1 Cell) The i-th entry of this cell
%                   array is a string which is the tag of the desired
%                   preset to execute as the i-th preset in the group. See
%                   main.m for a list of presets.
% group_settings    (Struct) This is where any group-shared settings (e.g.,
%                   system, initial conditions, etc.) are stored for
%                   subsequent initialization of each of the presets.
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

switch preset_group
    
    % *********************************************************************
    %
    % CASE STUDY I: 2ND ORDER ACADEMIC EXAMPLE, MINIMAL BASES
    %
    
    case 'CS_2ndorder_N1_3'
   
        %%   
 
        % *****************************************************************
        %
        % METHODS TO PERFORM THE SWEEP FOR
        %           
        
        alg_list =  {
                        'irl'
%                         'spi'
                        'radp_matched'
%                         'vi'
                                            };
 
        % Number of algorithms tested                
        numalgs = size(alg_list, 1);
       
        
        
        % *****************************************************************
        %
        % IC VECTORS TO SWEEP
        %   
        
        % Vector of x_1(0)
        x1vec = 1;
%         x1vec = 0.75;
%         x1vec = (-1:0.25:1)';
%         x1vec = (-1:0.5:1)';
%         x1vec = [-1;1];
%         x1vec = [-1;0;1];
%         x1vec = [0];
%         x1vec = [-1;0.5;1];        
        
        % Vector of x_2(0)
        x2vec = 1;
%         x2vec = -1;
%         x2vec = 0.5;
%         x2vec = (-1:0.25:1)';
%         x2vec = (-1:0.5:1)';
%         x2vec = [-1;1];
%         x2vec = [-1;0;1];
%         x2vec = [0];
%         x2vec = [-1;0.75;1];   

        % Size of IC vectors
        nx1 = size(x1vec,1);
        nx2 = size(x2vec,1);
        
        % Control to take out the IC (0,0) from the sweep, if desired
        remove_0_IC = 1;
        
        % Indices in the state vector corresponding to the state variables
        % swept
        xinds = [1; 2];
        
        % Store these IC vectors in a cell array for indexing
        x0_cell = [];
        ind_0_IC = -1;  % Index of the (0,0) IC (initialized to dummy val.)
        for i = 1:nx2
            x2 = x2vec(i);
            for j = 1:nx1
                x1 = x1vec(j);
                add_IC = ~(remove_0_IC && x1 == 0 && x2 == 0);
                if add_IC
                    x0tmp = zeros(2, 1);
                    x0tmp(xinds) = [x1 ; x2];
                    x0_cell = [ x0_cell ; {x0tmp} ];
                else
                    % Store index of the (0,0) IC
                    ind_0_IC = (i-1) * nx1 + j;
                end
            end
        end
        
        % Number of ICs tested in the sweep (should be nx1*nx2, or nx1*nx2
        % - 1 if 0 was taken out)
        numICs = size(x0_cell, 1);
        
        % Matrix of ICs to plot auxiliary plots for. E.g.,  if
        % aux_plot_x0s(2,:) = [1 -1], then the second set of ICs for which
        % the auxiliary plots will be plotted for will be for the IC x_0 =
        % [1; -1]
        aux_plot_x0s = [
%                         -1    -1 
%                         -1    1
                        1       -1
                        1       1
                        0.25    -0.5
                                ];
                            
        % *****************************************************************
        %
        % PRESET LIST
        %      
        
        % Total number of presets tested. Is the number of algorithms
        % tested times the number of ICs in the sweep; i.e., numalgs *
        % numICs
        numpresets = numalgs * numICs;
        
        preset_list = cell(numpresets, 1);
        
        % Initialize each preset in the list
        for i = 1:numalgs
            alg = alg_list{i};
            preset = [preset_group '_' alg];
            for j = 1:numICs 
                preset_list{numICs * (i - 1) + j} = preset; 
            end
        end

        % *****************************************************************
        %
        % PRESET GROUP SHARED SETTINGS      
        %        

        % System tag
        sys.tag = 'vrabie_lewis_2009_hard';
        
        % State penalty function
        Q = 'vrabie_lewis_2009_hard_ex';
        
        % Control penalty matrix
        R = 1;
        
        % Probing noise signal
        % See eval_noise.m for documentation of how the 'sum_sinusoids'
        % exploration noise option works.
        noise.tag = 'sum_sinusoids';
        
        % A cos(t)
        noise.cos1_sin0 = 1;
        noise.scalevec = 5;
%         noise.scalevec = 10;
%         noise.scalevec = 20;
        noise.freqvec = 1;
        
%         % A sin(t)
%         noise.cos1_sin0 = 0;
%         noise.scalevec = 5;
%         noise.freqvec = 1;

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 0; 1];
%         noise.scalevec = [5; 1; 1];
%         noise.freqvec = [1; 0.5; 0.1];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1];
%         noise.scalevec = [20; 1];
%         noise.freqvec = [1; 0.1];

        
        % Learning interval [0, t_f]
        tf = 10;
        
        % Time to simulate after learning [t_f, t_f + tsim]
        tsim = 0.001;
        
        % Number of iterations to terminate after (IRL and RADP only)
        istar = 5;
        
        % Initial stabilizing policy (if algorithm uses it)
%         u_0.tag = 'comp_vrabie_lewis_2009_hard';   
        u_0.tag = 'vrabie_lewis_2009_hard_min_critic'; 

        % ***********************
        %
        % CRITIC BASIS   
        %  
        
        % Default basis for critic NN; i.e., \hat{V}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        
%         basis_critic.tag = 'order_2_degree_4';        
        basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_debug';

        % ICs
        basis_critic.c0 = [0; 3/2; 3];
        
        % ***********************
        %
        % ACTOR BASIS -- g(x) UNKNOWN
        %  
        % Basis for actor NN without g(x); i.e., \hat{\mu}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        %
        
        basis_actor_no_g.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';   



        % ***********************
        %
        % HAMILTONIAN NETWORK BASIS (VI ONLY)
        %        
        
        % Hamiltonian NN basis \Theta(x)
%         Theta.tag = 'comp_vrabie_lewis_2009_hard_vi_theta';
        Theta.tag = 'comp_vrabie_lewis_2009_hard_vi_theta_w_Q';
        
        % Store parameters
        basis_hamiltonian.Theta = Theta;

        % ***********************
        %
        % OPTIMAL BASES, WEIGHTS
        %        
        
        % Tags corresponding to optimal bases, associated optimal weights
        basis_critic_opt.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
        basis_critic_opt.c_star = [0.5; 1; 1];
        
        basis_actor_no_g_opt.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';
        basis_actor_no_g_opt.w_star = [-1; -2];
        basis_actor_no_g_opt.w_star_VI = -2 * [-1; -2];
        
        basis_hamiltonian_opt.tag = 'comp_vrabie_lewis_2009_hard_vi_theta';
        basis_hamiltonian_opt.v_star = [1; 4; 4;];        
        
        
        % Optimal weights
        c_star = [0.5; 1; 1];
        w_star = [-1; -2]; 
        w_star_VI = -2 * w_star; 
%         v_star = [1; 4; 4];  
%         v_star = [1; 4; 4; -1]; 
        v_star = [1; 4; 4; -1; -1; -2];
        
        % Store if each of the bases used are indeed the optimal bases
        % (or contain the optimal basis functions as a subset)
        isopt_critic = 1;
        isopt_actor_no_g = 1;
        isopt_hamiltonian = 1;
        
        
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        
        
        % Do plots for each individual preset in the group
        do_individual_plots = 1;
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 1;
                      
        % Whether or not to plot the optimal weight values on surface plots
        % for critic weights
        plot_opt_weight = 1;
        
        
        % ***********************
        %
        % STORE SETTINGS
        %
        
        group_settings.Q = Q;
        group_settings.R = R;
        group_settings.noise = noise;
        group_settings.tf = tf;
        group_settings.tsim = tsim;
        
        group_settings.istar = istar;

        group_settings.alg_list = alg_list;
        group_settings.numalgs = numalgs;
        group_settings.x0_cell = x0_cell;
        group_settings.numICs = numICs;
        group_settings.x1vec = x1vec;
        group_settings.nx1 = nx1;
        group_settings.x2vec = x2vec;
        group_settings.nx2 = nx2;
        group_settings.xinds = xinds;
        group_settings.remove_0_IC = remove_0_IC;
        group_settings.ind_0_IC = ind_0_IC;
        
        group_settings.plot_opt_weight = plot_opt_weight;
        group_settings.c_star = c_star;
        group_settings.w_star = w_star;
        group_settings.w_star_VI = w_star_VI;
        group_settings.v_star = v_star;
        
        group_settings.basis_critic_opt = basis_critic_opt;
        group_settings.basis_actor_no_g_opt = basis_actor_no_g_opt;
        group_settings.basis_hamiltonian_opt = basis_hamiltonian_opt;
        
        group_settings.isopt_critic = isopt_critic;
        group_settings.isopt_actor_no_g = isopt_actor_no_g;
        group_settings.isopt_hamiltonian = isopt_hamiltonian;
        
        group_settings.aux_plot_x0s = aux_plot_x0s;

    % *********************************************************************
    %
    % CASE STUDY II: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 4
    % TERMS
    %
    
    case 'CS_2ndorder_N1_4'
   
        %%   
 
        % *****************************************************************
        %
        % METHODS TO PERFORM THE SWEEP FOR
        %           
        
        alg_list =  {
                        'irl'
                        'spi'
                        'radp_matched'
%                         'vi'
                                            };
 
        % Number of algorithms tested                
        numalgs = size(alg_list, 1);
       
        
        
        % *****************************************************************
        %
        % IC VECTORS TO SWEEP
        %   
        
        % Vector of x_1(0)
%         x1vec = 1;
%         x1vec = -1;
%         x1vec = 0.5;
        x1vec = (-1:0.25:1)';
%         x1vec = (-1:0.5:1)';
%         x1vec = [-1;1];
%         x1vec = [-1;0;1];
%         x1vec = [0];
%         x1vec = [-1;0.5;1];    
%         x1vec = [-1;-0.5;1];
        
        % Vector of x_2(0)
%         x2vec = 1;
%         x2vec = -1;
%         x2vec = -0.5;
%         x2vec = -0.75;
        x2vec = (-1:0.25:1)';
%         x2vec = (-1:0.5:1)';
%         x2vec = [-1;1];
%         x2vec = [-1;0;1];
%         x2vec = [0];
%         x2vec = [-1;-0.5;1];   

        % Size of IC vectors
        nx1 = size(x1vec,1);
        nx2 = size(x2vec,1);
        
        % Control to take out the IC (0,0) from the sweep, if desired
        remove_0_IC = 1;
        
        % Indices in the state vector corresponding to the state variables
        % swept
        xinds = [1; 2];
        
        % Store these IC vectors in a cell array for indexing
        x0_cell = [];
        ind_0_IC = -1;  % Index of the (0,0) IC (initialized to dummy val.)
        for i = 1:nx2
            x2 = x2vec(i);
            for j = 1:nx1
                x1 = x1vec(j);
                add_IC = ~(remove_0_IC && x1 == 0 && x2 == 0);
                if add_IC
                    x0tmp = zeros(2, 1);
                    x0tmp(xinds) = [x1 ; x2];
                    x0_cell = [ x0_cell ; {x0tmp} ];
                else
                    % Store index of the (0,0) IC
                    ind_0_IC = (i-1) * nx1 + j;
                end
            end
        end
        
        % Number of ICs tested in the sweep (should be nx1*nx2, or nx1*nx2
        % - 1 if 0 was taken out)
        numICs = size(x0_cell, 1);
        
        % Matrix of ICs to plot auxiliary plots for. E.g.,  if
        % aux_plot_x0s(2,:) = [1 -1], then the second set of ICs for which
        % the auxiliary plots will be plotted for will be for the IC x_0 =
        % [1; -1]
        aux_plot_x0s = [
%                         -1    -1 
%                         -1    1
                        1       -1
                        1       1
                        0.25    -0.5
                                ];
                            
        
        % *****************************************************************
        %
        % PRESET LIST
        %      
        
        % Total number of presets tested. Is the number of algorithms
        % tested times the number of ICs in the sweep; i.e., numalgs *
        % numICs
        numpresets = numalgs * numICs;
        
        preset_list = cell(numpresets, 1);
        
        % Initialize each preset in the list
        for i = 1:numalgs
            alg = alg_list{i};
            preset = [preset_group '_' alg];
            for j = 1:numICs 
                preset_list{numICs * (i - 1) + j} = preset; 
            end
        end

        % *****************************************************************
        %
        % PRESET GROUP SHARED SETTINGS      
        %        

        % System tag
        sys.tag = 'vrabie_lewis_2009_hard';
        
        % State penalty function
        Q = 'vrabie_lewis_2009_hard_ex';
        
        % Control penalty matrix
        R = 1;
        
        % Probing noise signal
        % See eval_noise.m for documentation of how the 'sum_sinusoids'
        % exploration noise option works.
        noise.tag = 'sum_sinusoids';
        
%         % A cos(t)
        noise.cos1_sin0 = 1;
%         noise.scalevec = 5;
%         noise.scalevec = 10;
        noise.scalevec = 20;
        noise.freqvec = 1;

%         % A cos(5t)
%         noise.cos1_sin0 = 1;
%         noise.scalevec = 5;
%         noise.freqvec = 5;
        
%         % A sin(t)
%         noise.cos1_sin0 = 0;
%         noise.scalevec = 5;
%         noise.freqvec = 1;

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1];
%         noise.scalevec = [20; -10];
%         noise.freqvec = [1; 0.1];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1];
%         noise.scalevec = [5; 1];
%         noise.freqvec = [1; 0.1];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 0; 1];
%         noise.scalevec = [5; 1; 1];
%         noise.freqvec = [1; 0.5; 0.1];
        
%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 0; 1];
%         noise.scalevec = [10; 5; 5];
%         noise.freqvec = [1; 0.5; 0.1];        

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1];
%         noise.scalevec = [5; 0.25];
%         noise.freqvec = [1; 0.1];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1; 1; 0];
%         noise.scalevec = [10; 1; 5; 5];
%         noise.freqvec = [1; 5; 0.1; 0.5];
        
        % Learning interval [0, t_f]
        tf = 10;
        
        % Time to simulate after learning [t_f, t_f + tsim]
        tsim = 10;
        
        % Number of iterations to terminate after (IRL and RADP only)
        istar = 5;
        
        % Initial stabilizing policy (if algorithm uses it)
%         u_0.tag = 'comp_vrabie_lewis_2009_hard';   
        u_0.tag = 'vrabie_lewis_2009_hard_min_critic'; 

        % ***********************
        %
        % CRITIC BASIS   
        %  
        
        % Default basis for critic NN; i.e., \hat{V}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        
        basis_critic.tag = 'vrabie_lewis_2009_hard_critic_nonmin';        
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_debug';
%         basis_critic.tag = 'order_2_degree_4';

        % ICs
        basis_critic.c0 = [0; 3/2; 3; 0];


        % ***********************
        %
        % ACTOR BASIS -- g(x) UNKNOWN
        %  
        % Basis for actor NN without g(x); i.e., \hat{\mu}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        %
        
        basis_actor_no_g.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';   


        % ***********************
        %
        % HAMILTONIAN NETWORK BASIS (VI ONLY)
        %        
        
        % Hamiltonian NN basis \Theta(x)
%         Theta.tag = 'comp_vrabie_lewis_2009_hard_vi_theta';
        Theta.tag = 'comp_vrabie_lewis_2009_hard_vi_theta_w_Q';
        
        % Store parameters
        basis_hamiltonian.Theta = Theta;

        % ***********************
        %
        % OPTIMAL BASES, WEIGHTS
        %        
        
        % Tags corresponding to optimal bases, associated optimal weights
        basis_critic_opt.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
        basis_critic_opt.c_star = [0.5; 1; 1];
        
        basis_actor_no_g_opt.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';
        basis_actor_no_g_opt.w_star = [-1; -2];
        basis_actor_no_g_opt.w_star_VI = -2 * [-1; -2];
        
        basis_hamiltonian_opt.tag = 'comp_vrabie_lewis_2009_hard_vi_theta';
        basis_hamiltonian_opt.v_star = [1; 4; 4;];
        
        % Optimal weights associated with the bases used (NOT with the
        % optimal bases)
        c_star = [0.5; 1; 1; 0];
        w_star = [-1; -2]; 
        w_star_VI = -2 * w_star; 
%         v_star = [1; 4; 4];  
%         v_star = [1; 4; 4; -1]; 
        v_star = [1; 4; 4; -1; -1; -2];  
        
        % Store if each of the bases used are indeed the optimal bases
        % (or contain the optimal basis functions as a subset)
        isopt_critic = 1;
        isopt_actor_no_g = 1;
        isopt_hamiltonian = 1;
        
        
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        
        
        % Do plots for each individual preset in the group
        do_individual_plots = 0;
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 1;
                      
        % Whether or not to plot the optimal weight values on surface plots
        % for critic weights
        plot_opt_weight = 1;

        
        % ***********************
        %
        % STORE SETTINGS
        %
        
        group_settings.Q = Q;
        group_settings.R = R;
        group_settings.noise = noise;
        group_settings.tf = tf;
        group_settings.tsim = tsim;
        
        group_settings.istar = istar;

        group_settings.alg_list = alg_list;
        group_settings.numalgs = numalgs;
        group_settings.x0_cell = x0_cell;
        group_settings.numICs = numICs;
        group_settings.x1vec = x1vec;
        group_settings.nx1 = nx1;
        group_settings.x2vec = x2vec;
        group_settings.nx2 = nx2;
        group_settings.xinds = xinds;
        group_settings.remove_0_IC = remove_0_IC;
        group_settings.ind_0_IC = ind_0_IC;
        
        group_settings.plot_opt_weight = plot_opt_weight;
        group_settings.c_star = c_star;
        group_settings.w_star = w_star;
        group_settings.w_star_VI = w_star_VI;
        group_settings.v_star = v_star;

        group_settings.basis_critic_opt = basis_critic_opt;
        group_settings.basis_actor_no_g_opt = basis_actor_no_g_opt;
        group_settings.basis_hamiltonian_opt = basis_hamiltonian_opt;
        
        group_settings.isopt_critic = isopt_critic;
        group_settings.isopt_actor_no_g = isopt_actor_no_g;
        group_settings.isopt_hamiltonian = isopt_hamiltonian;
        
        group_settings.aux_plot_x0s = aux_plot_x0s;
            
 
       
        
    % *********************************************************************
    %
    % CASE STUDY III: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_2^4 REMOVED
    %
    
    case 'CS_2ndorder_N1_7_no_x2p4'
   
        %%   
 
        % *****************************************************************
        %
        % METHODS TO PERFORM THE SWEEP FOR
        %           
        
        alg_list =  {
%                         'irl'
                        'spi'
%                         'radp_matched'
%                         'vi'
                                            };
 
        % Number of algorithms tested                
        numalgs = size(alg_list, 1);
       
        
        
        % *****************************************************************
        %
        % IC VECTORS TO SWEEP
        %   
        
        % Vector of x_1(0)
        x1vec = 1;
%         x1vec = -1;
%         x1vec = 0.25;
%         x1vec = 0.01;
%         x1vec = (-1:0.25:1)';
%         x1vec = (-1:0.5:1)';
%         x1vec = [-1;1];
%         x1vec = [0.99;1];
%         x1vec = [-1;0;1];
%         x1vec = [0];
%         x1vec = [-1;0.5;1];        
        
        % Vector of x_2(0)
        x2vec = 1;
%         x2vec = -1;
%         x2vec = 0.25;
%         x2vec = 0.01;
%         x2vec = -0.5;
%         x2vec = 0.75;
%         x2vec = (-1:0.25:1)';
%         x2vec = (-1:0.5:1)';
%         x2vec = [-1;1];
%         x2vec = [0.99;1];
%         x2vec = [-1;0;1];
%         x2vec = [0];
%         x2vec = [-1;0.75;1];   

        % Size of IC vectors
        nx1 = size(x1vec,1);
        nx2 = size(x2vec,1);
        
        % Control to take out the IC (0,0) from the sweep, if desired
        remove_0_IC = 1;
        
        % Indices in the state vector corresponding to the state variables
        % swept
        xinds = [1; 2];
        
        % Store these IC vectors in a cell array for indexing
        x0_cell = [];
        ind_0_IC = -1;  % Index of the (0,0) IC (initialized to dummy val.)
        for i = 1:nx2
            x2 = x2vec(i);
            for j = 1:nx1
                x1 = x1vec(j);
                add_IC = ~(remove_0_IC && x1 == 0 && x2 == 0);
                if add_IC
                    x0tmp = zeros(2, 1);
                    x0tmp(xinds) = [x1 ; x2];
                    x0_cell = [ x0_cell ; {x0tmp} ];
                else
                    % Store index of the (0,0) IC
                    ind_0_IC = (i-1) * nx1 + j;
                end
            end
        end
        
        % Number of ICs tested in the sweep (should be nx1*nx2, or nx1*nx2
        % - 1 if 0 was taken out)
        numICs = size(x0_cell, 1);
        
        % Matrix of ICs to plot auxiliary plots for. E.g.,  if
        % aux_plot_x0s(2,:) = [1 -1], then the second set of ICs for which
        % the auxiliary plots will be plotted for will be for the IC x_0 =
        % [1; -1]
        aux_plot_x0s = [
%                         -1    -1 
%                         -1    1
                        1       -1
                        1       1
                        0.25    -0.5
%                         0.25    0.75
                                ];
                            
               
        
        % *****************************************************************
        %
        % PRESET LIST
        %      
        
        % Total number of presets tested. Is the number of algorithms
        % tested times the number of ICs in the sweep; i.e., numalgs *
        % numICs
        numpresets = numalgs * numICs;
        
        preset_list = cell(numpresets, 1);
        
        % Initialize each preset in the list
        for i = 1:numalgs
            alg = alg_list{i};
            preset = [preset_group '_' alg];
            for j = 1:numICs 
                preset_list{numICs * (i - 1) + j} = preset; 
            end
        end

        % *****************************************************************
        %
        % PRESET GROUP SHARED SETTINGS      
        %        

        % System tag
        sys.tag = 'vrabie_lewis_2009_hard';
        
        % State penalty function
        Q = 'vrabie_lewis_2009_hard_ex';
        
        % Control penalty matrix
        R = 1;
        
        % Probing noise signal
        % See eval_noise.m for documentation of how the 'sum_sinusoids'
        % exploration noise option works.
        
%         % Decaying sinusoid
%         noise.tag = 'dacaying_sinusoid';
        
        % Sum of sinusoids
        noise.tag = 'sum_sinusoids';
        
        % A cos(t)
        noise.cos1_sin0 = 1;
        noise.scalevec = 3.5;
%         noise.scalevec = 5;
%         noise.scalevec = 10;
%         noise.scalevec = 20;
        noise.freqvec = 1;

%         % A cos(5t)
%         noise.cos1_sin0 = 1;
%         noise.scalevec = 5;
%         noise.freqvec = 5;
        
%         % A sin(t)
%         noise.cos1_sin0 = 0;
%         noise.scalevec = 5;
%         noise.freqvec = 1;

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1];
%         noise.scalevec = [5; 1];
%         noise.freqvec = [1; 0.1];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 0; 1];
%         noise.scalevec = [5; 1; 1];
%         noise.freqvec = [1; 0.5; 0.1];
        
%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 0; 1];
%         noise.scalevec = [10; 5; 5];
%         noise.freqvec = [1; 0.5; 0.1];        

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1; 0 ];
%         noise.scalevec = [20; 5; 5 ];
%         noise.freqvec = [1; 5; 0.1 ];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1; 1; 0];
%         noise.scalevec = [10; 1; 5; 5];
%         noise.freqvec = [1; 5; 0.1; 0.5];
        
        % Learning interval [0, t_f]
        tf = 10;
        
        % Time to simulate after learning [t_f, t_f + tsim]
        tsim = 10;
        
        % Number of iterations to terminate after (IRL and RADP only)
        istar = 15;
        
        % Initial stabilizing policy (if algorithm uses it)
%         u_0.tag = 'comp_vrabie_lewis_2009_hard';   
%         u_0.tag = 'vrabie_lewis_2009_hard_min_critic'; 
        u_0.tag = 'vrabie_lewis_2009_hard_x_2_only'; 
%         u_0.tag = 'vrabie_lewis_2009_hard_testing'; 

        % ***********************
        %
        % CRITIC BASIS   
        %  
        
        % Default basis for critic NN; i.e., \hat{V}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        
        basis_critic.tag = 'vrabie_lewis_2009_hard_critic_7terms_no_x_2^4';
%         basis_critic.tag = 'vrabie_lewis_2009_hard_critic_8terms';
%         basis_critic.tag = 'vrabie_lewis_2009_hard_critic_nonmin';         
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_debug';
%         basis_critic.tag = 'order_2_degree_4';


        % ICs
%         basis_critic.c0 = [0 ; 5/2; zeros(5,1)];
%         basis_critic.c0 = [0 ; 10/2; 0; 0; 0; -0.2/2; 0];
        basis_critic.c0 = [0 ; 10/2; 0; 0; 0; 0; 0];
        
%         basis_critic.c0 = [0.1871
%         1.8123
%         0.0000
%         0.2025
%         -0.0000
%         0.1876
%         -0.0000];

%         basis_critic.c0 = [  0.7708
%                                     1.9715
%                                     0.3448
%                                     -0.2336
%                                     -0.1538
%                                     -0.0826
%                                     -0.1790  ];


%         basis_critic.c0 = [ 0.6821
%             1.8514
%             0.2875
%            -0.2727
%            -0.1916
%            -0.1255
%            -0.2007 ];

        % ***********************
        %
        % ACTOR BASIS -- g(x) UNKNOWN
        %  
        % Basis for actor NN without g(x); i.e., \hat{\mu}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        %
        
        basis_actor_no_g.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';   



        % ***********************
        %
        % HAMILTONIAN NETWORK BASIS (VI ONLY)
        %        
        
        % Hamiltonian NN basis \Theta(x)
%         Theta.tag = 'comp_vrabie_lewis_2009_hard_vi_theta';
        Theta.tag = 'comp_vrabie_lewis_2009_hard_vi_theta_w_Q';
        
        % Store parameters
        basis_hamiltonian.Theta = Theta;

        % ***********************
        %
        % OPTIMAL BASES, WEIGHTS
        %        
        
        % Tags corresponding to optimal bases, associated optimal weights
        basis_critic_opt.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
        basis_critic_opt.c_star = [0.5; 1; 1];
        
        basis_actor_no_g_opt.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';
        basis_actor_no_g_opt.w_star = [-1; -2];
        basis_actor_no_g_opt.w_star_VI = -2 * [-1; -2];
        
        basis_hamiltonian_opt.tag = 'comp_vrabie_lewis_2009_hard_vi_theta';
        basis_hamiltonian_opt.v_star = [1; 4; 4;];
        
        % Optimal weights associated with the bases used (NOT with the
        % optimal bases)
%         c_star = [0.5; 1; 1; 0];
        w_star = [-1; -2]; 
        w_star_VI = -2 * w_star; 
%         v_star = [1; 4; 4];  
%         v_star = [1; 4; 4; -1]; 
        v_star = [1; 4; 4; -1; -1; -2];  
        
        % Store if each of the bases used are indeed the optimal bases
        % (or contain the optimal basis functions as a subset)
        isopt_critic = 0;
        isopt_actor_no_g = 1;
        isopt_hamiltonian = 1;
        
        
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        
        
        % Do plots for each individual preset in the group
        do_individual_plots = 0;
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 1;
                      
        % Whether or not to plot the optimal weight values on surface plots
        % for critic weights
        plot_opt_weight = 1;

        
        % ***********************
        %
        % STORE SETTINGS
        %
        
        group_settings.Q = Q;
        group_settings.R = R;
        group_settings.noise = noise;
        group_settings.tf = tf;
        group_settings.tsim = tsim;
        
        group_settings.istar = istar;

        group_settings.alg_list = alg_list;
        group_settings.numalgs = numalgs;
        group_settings.x0_cell = x0_cell;
        group_settings.numICs = numICs;
        group_settings.x1vec = x1vec;
        group_settings.nx1 = nx1;
        group_settings.x2vec = x2vec;
        group_settings.nx2 = nx2;
        group_settings.xinds = xinds;
        group_settings.remove_0_IC = remove_0_IC;
        group_settings.ind_0_IC = ind_0_IC;
        
        group_settings.plot_opt_weight = plot_opt_weight;
%         group_settings.c_star = c_star;
        group_settings.w_star = w_star;
        group_settings.w_star_VI = w_star_VI;
        group_settings.v_star = v_star;

        group_settings.basis_critic_opt = basis_critic_opt;
        group_settings.basis_actor_no_g_opt = basis_actor_no_g_opt;
        group_settings.basis_hamiltonian_opt = basis_hamiltonian_opt;
        
        group_settings.isopt_critic = isopt_critic;
        group_settings.isopt_actor_no_g = isopt_actor_no_g;
        group_settings.isopt_hamiltonian = isopt_hamiltonian;
        
        group_settings.aux_plot_x0s = aux_plot_x0s;             

        
    % *********************************************************************
    %
    % CASE STUDY IV: CART INVERTED PENDULUM SYSTEM
    %
    
    case 'CS_cip'
   
        %%   
 
        % *****************************************************************
        %
        % METHODS TO PERFORM THE SWEEP FOR
        %           
        
        alg_list =  {
%                         'irl'
%                         'spi'
%                         'radp_matched'
                        'vi'
                                            };
 
        % Number of algorithms tested                
        numalgs = size(alg_list, 1);
       
        
        
        % *****************************************************************
        %
        % IC VECTORS TO SWEEP
        %   
        
        % Vector of x_1(0)
        x1vec = 1;
%         x1vec = -1;
%         x1vec = 0.25;
%         x1vec = 0.01;
%         x1vec = (-1:0.25:1)';
%         x1vec = (-1:0.5:1)';
%         x1vec = [-1;1];
%         x1vec = [0.99;1];
%         x1vec = [-1;0;1];
%         x1vec = [0];
%         x1vec = [-1;0.5;1];        
        
        % Vector of x_2(0)
        x2vec = 15 * pi/180;
%         x2vec = 1;
%         x2vec = -1;
%         x2vec = 0.25;
%         x2vec = 0.01;
%         x2vec = 0.5;
%         x2vec = 0.75;
%         x2vec = 0.5 * (-1:0.25:1)';
%         x2vec = 0.5 * (-1:0.5:1)';
%         x2vec = [-1;1];
%         x2vec = [0.99;1];
%         x2vec = 0.5 * [-1;0;1];
%         x2vec = [0];
%         x2vec = [-1;0.75;1];   

        % Size of IC vectors
        nx1 = size(x1vec,1);
        nx2 = size(x2vec,1);
        
        % Control to take out the IC (0,0) from the sweep, if desired
        remove_0_IC = 1;
        
        % Indices in the state vector corresponding to the state variables
        % swept
        xinds = [1; 3];
        
        % Store these IC vectors in a cell array for indexing
        x0_cell = [];
        ind_0_IC = -1;  % Index of the (0,0) IC (initialized to dummy val.)
        for i = 1:nx2
            x2 = x2vec(i);
            for j = 1:nx1
                x1 = x1vec(j);
                add_IC = ~(remove_0_IC && x1 == 0 && x2 == 0);
                if add_IC
                    x0tmp = zeros(4, 1);
                    x0tmp(xinds) = [x1 ; x2];
                    x0_cell = [ x0_cell ; {x0tmp} ];
                else
                    % Store index of the (0,0) IC
                    ind_0_IC = (i-1) * nx1 + j;
                end
            end
        end
        
        % Number of ICs tested in the sweep (should be nx1*nx2, or nx1*nx2
        % - 1 if 0 was taken out)
        numICs = size(x0_cell, 1);
        
        % Matrix of ICs to plot auxiliary plots for. E.g.,  if
        % aux_plot_x0s(2,:) = [1 -1], then the second set of ICs for which
        % the auxiliary plots will be plotted for will be for the IC x_0 =
        % [1; -1]
%         aux_plot_x0s = [
% %                         -1      0   -1      0   
% %                         -1      0   1       0
%                         1       0   0.5 	0
% %                         1       0   1       0      
% %                         0.25    0   0.25    0
% %                         0.25    0   0.75    0
%                                 ];
        aux_plot_x0s = zeros(1,4);
        aux_plot_x0s(xinds) = [x1vec ; x2vec];
                            
               
        
        % *****************************************************************
        %
        % PRESET LIST
        %      
        
        % Total number of presets tested. Is the number of algorithms
        % tested times the number of ICs in the sweep; i.e., numalgs *
        % numICs
        numpresets = numalgs * numICs;
        
        preset_list = cell(numpresets, 1);
        
        % Initialize each preset in the list
        for i = 1:numalgs
            alg = alg_list{i};
            preset = [preset_group '_' alg];
            for j = 1:numICs 
                preset_list{numICs * (i - 1) + j} = preset; 
            end
        end

        % *****************************************************************
        %
        % PRESET GROUP SHARED SETTINGS      
        %        

        % System tag
        sys.tag = 'cip';
        
        % State penalty function
        Q = eye(4);
        
        % Control penalty matrix
        R = 1;
        
        % Probing noise signal
        % See eval_noise.m for documentation of how the 'sum_sinusoids'
        % exploration noise option works.
        
%         % Decaying sinusoid
%         noise.tag = 'dacaying_sinusoid';
        
        % Sum of sinusoids
        noise.tag = 'sum_sinusoids';
        
        % A cos(t)
        noise.cos1_sin0 = 1;
%         noise.scalevec = 7.5;
        noise.scalevec = 5;
%         noise.scalevec = 10;
%         noise.scalevec = 20;
        noise.freqvec = 1;

%         % A cos(5t)
%         noise.cos1_sin0 = 1;
%         noise.scalevec = 5;
%         noise.freqvec = 5;
        
%         % A sin(t)
%         noise.cos1_sin0 = 0;
%         noise.scalevec = 5;
%         noise.freqvec = 1;

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1];
%         noise.scalevec = [5; 1];
%         noise.freqvec = [1; 0.1];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 0; 1];
%         noise.scalevec = [5; 1; 1];
%         noise.freqvec = [1; 0.5; 0.1];
        
%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 0; 1];
%         noise.scalevec = [10; 5; 5];
%         noise.freqvec = [1; 0.5; 0.1];        

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1; 0 ];
%         noise.scalevec = [20; 5; 5 ];
%         noise.freqvec = [1; 5; 0.1 ];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1; 1; 0];
%         noise.scalevec = [10; 1; 5; 5];
%         noise.freqvec = [1; 5; 0.1; 0.5];
        
        % Learning interval [0, t_f]
        tf = 10;
        
        % Time to simulate after learning [t_f, t_f + tsim]
        tsim = 0.01;
        
        % Number of iterations to terminate after (IRL and RADP only)
        istar = 20;
        
        % Initial stabilizing policy (if algorithm uses it)
        u_0.tag = 'actor_g_known_lqr';
%         u_0.tag = 'comp_vrabie_lewis_2009_hard';   
%         u_0.tag = 'vrabie_lewis_2009_hard_min_critic'; 
%         u_0.tag = 'vrabie_lewis_2009_hard_x_2_only'; 
%         u_0.tag = 'vrabie_lewis_2009_hard_testing'; 

        % ***********************
        %
        % CRITIC BASIS   
        %  
        
        % Default basis for critic NN; i.e., \hat{V}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
       
        % All monomials of total degree = 2
        basis_critic.tag = 'monomial';
        basis_critic.type = 'total_deg_K';   
        basis_critic.K = 2;
        basis_critic.IC_type = 'lqr';
        
%         % All even monomials of total degree <= 4
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'even';   
%         basis_critic.K = 4;
%         basis_critic.IC_type = 'lqr';

%         basis_critic.tag = 'vrabie_lewis_2009_hard_critic_7terms_no_x_2^4';
%         basis_critic.tag = 'vrabie_lewis_2009_hard_critic_8terms';
%         basis_critic.tag = 'vrabie_lewis_2009_hard_critic_nonmin';         
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_debug';
%         basis_critic.tag = 'order_2_degree_4';        
        

        % ***********************
        %
        % ACTOR BASIS -- g(x) UNKNOWN
        %  
        % Basis for actor NN without g(x); i.e., \hat{\mu}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        %
        
        % All odd monomials of total degree <= 3
        basis_actor_no_g.tag = 'monomial';
        basis_actor_no_g.type = 'odd';   
        basis_actor_no_g.K = 3;
        basis_actor_no_g.IC_type = 'none';

        % ***********************
        %
        % HAMILTONIAN NETWORK BASIS (VI ONLY)
        %        
        
        % Hamiltonian NN basis \Theta(x)
        
        % All monomials of total degree = 2
        Theta.tag = 'monomial';
        Theta.type = 'total_deg_K';   
        Theta.K = 2;
        Theta.IC_type = 'lqr';
        
        % Store parameters
        basis_hamiltonian.Theta = Theta;

        % ***********************
        %
        % OPTIMAL BASES, WEIGHTS
        %        
        
        
        % Store if each of the bases used are indeed the optimal bases
        % (or contain the optimal basis functions as a subset)
        isopt_critic = 0;
        isopt_actor_no_g = 0;
        isopt_hamiltonian = 0;        
        
        
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        
        
        % Do plots for each individual preset in the group
        do_individual_plots = 1;
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 1;
                      
        % Whether or not to plot the optimal weight values on surface plots
        % for critic weights
        plot_opt_weight = 0;

        
        % ***********************
        %
        % STORE SETTINGS
        %
        
        group_settings.Q = Q;
        group_settings.R = R;
        group_settings.noise = noise;
        group_settings.tf = tf;
        group_settings.tsim = tsim;
        
        group_settings.istar = istar;

        group_settings.alg_list = alg_list;
        group_settings.numalgs = numalgs;
        group_settings.x0_cell = x0_cell;
        group_settings.numICs = numICs;
        group_settings.x1vec = x1vec;
        group_settings.nx1 = nx1;
        group_settings.x2vec = x2vec;
        group_settings.nx2 = nx2;
        group_settings.xinds = xinds;
        group_settings.remove_0_IC = remove_0_IC;
        group_settings.ind_0_IC = ind_0_IC;
        
        group_settings.plot_opt_weight = plot_opt_weight;

        group_settings.isopt_critic = isopt_critic;
        group_settings.isopt_actor_no_g = isopt_actor_no_g;
        group_settings.isopt_hamiltonian = isopt_hamiltonian;        
        
        group_settings.aux_plot_x0s = aux_plot_x0s;             
        
        
        
        
        
    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 8
    % TERMS
    %
    
    case 'CS_2ndorder_N1_8'
   
        %%   
 
        % *****************************************************************
        %
        % METHODS TO PERFORM THE SWEEP FOR
        %           
        
        alg_list =  {
                        'irl'
                        'spi'
                        'radp_matched'
                        'vi'
                                            };
 
        % Number of algorithms tested                
        numalgs = size(alg_list, 1);
       
        
        
        % *****************************************************************
        %
        % IC VECTORS TO SWEEP
        %   
        
        % Vector of x_1(0)
%         x1vec = 1;
%         x1vec = -1;
%         x1vec = 0.75;
%         x1vec = (-1:0.25:1)';
%         x1vec = (-1:0.5:1)';
        x1vec = [-1;1];
%         x1vec = [-1;0;1];
%         x1vec = [0];
%         x1vec = [-1;-0.5;1];        
        
        % Vector of x_2(0)
%         x2vec = 1;
%         x2vec = -1;
%         x2vec = 0.5;
%         x2vec = 0.25;
%         x2vec = (-1:0.25:1)';
%         x2vec = (-1:0.5:1)';
        x2vec = [-1;1];
%         x2vec = [-1;0;1];
%         x2vec = [0];
%         x2vec = [-1;0.75;1];   

        % Size of IC vectors
        nx1 = size(x1vec,1);
        nx2 = size(x2vec,1);
        
        % Control to take out the IC (0,0) from the sweep, if desired
        remove_0_IC = 1;
        
        % Indices in the state vector corresponding to the state variables
        % swept
        xinds = [1; 2];
        
        % Store these IC vectors in a cell array for indexing
        x0_cell = [];
        ind_0_IC = -1;  % Index of the (0,0) IC (initialized to dummy val.)
        for i = 1:nx2
            x2 = x2vec(i);
            for j = 1:nx1
                x1 = x1vec(j);
                add_IC = ~(remove_0_IC && x1 == 0 && x2 == 0);
                if add_IC
                    x0tmp = zeros(2, 1);
                    x0tmp(xinds) = [x1 ; x2];
                    x0_cell = [ x0_cell ; {x0tmp} ];
                else
                    % Store index of the (0,0) IC
                    ind_0_IC = (i-1) * nx1 + j;
                end
            end
        end
        
        % Number of ICs tested in the sweep (should be nx1*nx2, or nx1*nx2
        % - 1 if 0 was taken out)
        numICs = size(x0_cell, 1);
        
        % Matrix of ICs to plot auxiliary plots for. E.g.,  if
        % aux_plot_x0s(2,:) = [1 -1], then the second set of ICs for which
        % the auxiliary plots will be plotted for will be for the IC x_0 =
        % [1; -1]
        aux_plot_x0s = [
                        -1      -1 
%                         -0.5    -1
                        -1      1
                        1       -1
                        1       1
                                    ];
                            
        
        % *****************************************************************
        %
        % PRESET LIST
        %      
        
        % Total number of presets tested. Is the number of algorithms
        % tested times the number of ICs in the sweep; i.e., numalgs *
        % numICs
        numpresets = numalgs * numICs;
        
        preset_list = cell(numpresets, 1);
        
        % Initialize each preset in the list
        for i = 1:numalgs
            alg = alg_list{i};
            preset = [preset_group '_' alg];
            for j = 1:numICs 
                preset_list{numICs * (i - 1) + j} = preset; 
            end
        end

        % *****************************************************************
        %
        % PRESET GROUP SHARED SETTINGS      
        %        

        % System tag
        sys.tag = 'vrabie_lewis_2009_hard';
        
        % State penalty function
        Q = 'vrabie_lewis_2009_hard_ex';
        
        % Control penalty matrix
        R = 1;
        
        % Probing noise signal
        % See eval_noise.m for documentation of how the 'sum_sinusoids'
        % exploration noise option works.
        noise.tag = 'sum_sinusoids';
        
%         % A cos(t)
%         noise.cos1_sin0 = 1;
%         noise.scalevec = 5;
% %         noise.scalevec = 10;
% %         noise.scalevec = 20;
%         noise.freqvec = 1;

%         % A cos(5t)
%         noise.cos1_sin0 = 1;
%         noise.scalevec = 5;
%         noise.freqvec = 5;
        
%         % A sin(t)
%         noise.cos1_sin0 = 0;
%         noise.scalevec = 5;
%         noise.freqvec = 1;

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1];
%         noise.scalevec = [5; 1];
%         noise.freqvec = [1; 0.1];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 0; 1];
%         noise.scalevec = [5; 1; 1];
%         noise.freqvec = [1; 0.5; 0.1];
        
        % Multi-sinusoid signal
        noise.cos1_sin0 = [1; 0; 1];
        noise.scalevec = [10; 5; 5];
        noise.freqvec = [1; 0.5; 0.1];        

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1; 0 ];
%         noise.scalevec = [20; 5; 5 ];
%         noise.freqvec = [1; 5; 0.1 ];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1; 1; 0];
%         noise.scalevec = [10; 1; 5; 5];
%         noise.freqvec = [1; 5; 0.1; 0.5];
        
        % Learning interval [0, t_f]
        tf = 5;
        
        % Time to simulate after learning [t_f, t_f + tsim]
        tsim = 10;
        
        % Number of iterations to terminate after (IRL and RADP only)
        istar = 5;
        
        % Initial stabilizing policy (if algorithm uses it)
%         u_0.tag = 'comp_vrabie_lewis_2009_hard';   
        u_0.tag = 'vrabie_lewis_2009_hard_min_critic'; 

        % ***********************
        %
        % CRITIC BASIS   
        %  
        
        % Default basis for critic NN; i.e., \hat{V}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        
        basis_critic.tag = 'vrabie_lewis_2009_hard_critic_8terms';
%         basis_critic.tag = 'vrabie_lewis_2009_hard_critic_nonmin';         
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_debug';
%         basis_critic.tag = 'order_2_degree_4';

%         % All even monomials of total degree <= 4
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'even';   
%         basis_critic.K = 4;
%         basis_critic.IC_type = 'none';

%         % All monomials of total degree = 2
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'total_deg_K';   
%         basis_critic.K = 2;
%         basis_critic.IC_type = 'none';

        % ICs
        basis_critic.c0 = [0; 3/2; 3; zeros(5,1)];
        
        % ***********************
        %
        % ACTOR BASIS -- g(x) UNKNOWN
        %  
        % Basis for actor NN without g(x); i.e., \hat{\mu}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        %
        
        basis_actor_no_g.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';   

%         % All even monomials of total degree <= 2
%         basis_actor_no_g.tag = 'monomial';
%         basis_actor_no_g.type = 'even';   
%         basis_actor_no_g.K = 2;
%         basis_actor_no_g.IC_type = 'none';

        % ***********************
        %
        % HAMILTONIAN NETWORK BASIS (VI ONLY)
        %        
        
        % Hamiltonian NN basis \Theta(x)
%         Theta.tag = 'comp_vrabie_lewis_2009_hard_vi_theta';
        Theta.tag = 'comp_vrabie_lewis_2009_hard_vi_theta_w_Q';
        
        % Store parameters
        basis_hamiltonian.Theta = Theta;

        % ***********************
        %
        % OPTIMAL BASES, WEIGHTS
        %        
        
        % Tags corresponding to optimal bases, associated optimal weights
        basis_critic_opt.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
        basis_critic_opt.c_star = [0.5; 1; 1];
        
        basis_actor_no_g_opt.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';
        basis_actor_no_g_opt.w_star = [-1; -2];
        basis_actor_no_g_opt.w_star_VI = -2 * [-1; -2];
        
        basis_hamiltonian_opt.tag = 'comp_vrabie_lewis_2009_hard_vi_theta';
        basis_hamiltonian_opt.v_star = [1; 4; 4;];
        
        % Optimal weights associated with the bases used (NOT with the
        % optimal bases)
        c_star = [0.5; 1; 1; 0];
        w_star = [-1; -2]; 
        w_star_VI = -2 * w_star; 
        v_star = [1; 4; 4; -1];  
        
        % Store if each of the bases used are indeed the optimal bases
        % (or contain the optimal basis functions as a subset)
        isopt_critic = 1;
        isopt_actor_no_g = 1;
        isopt_hamiltonian = 1;
        
        
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        
        
        % Do plots for each individual preset in the group
        do_individual_plots = 0;
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 1;
                      
        % Whether or not to plot the optimal weight values on surface plots
        % for critic weights
        plot_opt_weight = 1;

        
        % ***********************
        %
        % STORE SETTINGS
        %
        
        group_settings.Q = Q;
        group_settings.R = R;
        group_settings.noise = noise;
        group_settings.tf = tf;
        group_settings.tsim = tsim;
        
        group_settings.istar = istar;

        group_settings.alg_list = alg_list;
        group_settings.numalgs = numalgs;
        group_settings.x0_cell = x0_cell;
        group_settings.numICs = numICs;
        group_settings.x1vec = x1vec;
        group_settings.nx1 = nx1;
        group_settings.x2vec = x2vec;
        group_settings.nx2 = nx2;
        group_settings.xinds = xinds;
        group_settings.remove_0_IC = remove_0_IC;
        group_settings.ind_0_IC = ind_0_IC;
        
        group_settings.plot_opt_weight = plot_opt_weight;
        group_settings.c_star = c_star;
        group_settings.w_star = w_star;
        group_settings.w_star_VI = w_star_VI;
        group_settings.v_star = v_star;

        group_settings.basis_critic_opt = basis_critic_opt;
        group_settings.basis_actor_no_g_opt = basis_actor_no_g_opt;
        group_settings.basis_hamiltonian_opt = basis_hamiltonian_opt;
        
        group_settings.isopt_critic = isopt_critic;
        group_settings.isopt_actor_no_g = isopt_actor_no_g;
        group_settings.isopt_hamiltonian = isopt_hamiltonian;
        
        group_settings.aux_plot_x0s = aux_plot_x0s;             
        
        

    % *********************************************************************
    %
    % CASE STUDY: 2ND ORDER ACADEMIC EXAMPLE, CRITIC BASIS WITH N_1 = 7
    % TERMS, ESSENTIAL ACTIVATION FUNCTION x_1^2 REMOVED
    %
    
    case 'CS_2ndorder_N1_7_no_x1p2'
   
        %%   
 
        % *****************************************************************
        %
        % METHODS TO PERFORM THE SWEEP FOR
        %           
        
        alg_list =  {
                        'irl'
%                         'spi'
%                         'radp_matched'
%                         'vi'
                                            };
 
        % Number of algorithms tested                
        numalgs = size(alg_list, 1);
       
        
        
        % *****************************************************************
        %
        % IC VECTORS TO SWEEP
        %   
        
        % Vector of x_1(0)
        x1vec = 1;
%         x1vec = -1;
%         x1vec = 0.75;
%         x1vec = (-1:0.25:1)';
%         x1vec = (-1:0.5:1)';
%         x1vec = [-1;1];
%         x1vec = [-1;0;1];
%         x1vec = [0];
%         x1vec = [-1;0.5;1];        
        
        % Vector of x_2(0)
        x2vec = 1;
%         x2vec = -1;
%         x2vec = 0.5;
%         x2vec = 0.25;
%         x2vec = (-1:0.25:1)';
%         x2vec = (-1:0.5:1)';
%         x2vec = [-1;1];
%         x2vec = [-1;0;1];
%         x2vec = [0];
%         x2vec = [-1;0.75;1];   

        % Size of IC vectors
        nx1 = size(x1vec,1);
        nx2 = size(x2vec,1);
        
        % Control to take out the IC (0,0) from the sweep, if desired
        remove_0_IC = 1;
        
        % Indices in the state vector corresponding to the state variables
        % swept
        xinds = [1; 2];
        
        % Store these IC vectors in a cell array for indexing
        x0_cell = [];
        ind_0_IC = -1;  % Index of the (0,0) IC (initialized to dummy val.)
        for i = 1:nx2
            x2 = x2vec(i);
            for j = 1:nx1
                x1 = x1vec(j);
                add_IC = ~(remove_0_IC && x1 == 0 && x2 == 0);
                if add_IC
                    x0tmp = zeros(2, 1);
                    x0tmp(xinds) = [x1 ; x2];
                    x0_cell = [ x0_cell ; {x0tmp} ];
                else
                    % Store index of the (0,0) IC
                    ind_0_IC = (i-1) * nx1 + j;
                end
            end
        end
        
        % Number of ICs tested in the sweep (should be nx1*nx2, or nx1*nx2
        % - 1 if 0 was taken out)
        numICs = size(x0_cell, 1);
        
        % Matrix of ICs to plot auxiliary plots for. E.g.,  if
        % aux_plot_x0s(2,:) = [1 -1], then the second set of ICs for which
        % the auxiliary plots will be plotted for will be for the IC x_0 =
        % [1; -1]
        aux_plot_x0s = [
                        -1  -1 
                        -1  1
                        1   -1
                        1   1
                                ];
                            
        
        % *****************************************************************
        %
        % PRESET LIST
        %      
        
        % Total number of presets tested. Is the number of algorithms
        % tested times the number of ICs in the sweep; i.e., numalgs *
        % numICs
        numpresets = numalgs * numICs;
        
        preset_list = cell(numpresets, 1);
        
        % Initialize each preset in the list
        for i = 1:numalgs
            alg = alg_list{i};
            preset = [preset_group '_' alg];
            for j = 1:numICs 
                preset_list{numICs * (i - 1) + j} = preset; 
            end
        end

        % *****************************************************************
        %
        % PRESET GROUP SHARED SETTINGS      
        %        

        % System tag
        sys.tag = 'vrabie_lewis_2009_hard';
        
        % State penalty function
        Q = 'vrabie_lewis_2009_hard_ex';
        
        % Control penalty matrix
        R = 1;
        
        % Probing noise signal
        % See eval_noise.m for documentation of how the 'sum_sinusoids'
        % exploration noise option works.
        noise.tag = 'sum_sinusoids';
        
        % A cos(t)
        noise.cos1_sin0 = 1;
%         noise.scalevec = 3;
%         noise.scalevec = 5;
        noise.scalevec = 10;
%         noise.scalevec = 20;
        noise.freqvec = 1;

%         % A cos(5t)
%         noise.cos1_sin0 = 1;
%         noise.scalevec = 5;
%         noise.freqvec = 5;
        
%         % A sin(t)
%         noise.cos1_sin0 = 0;
%         noise.scalevec = 5;
%         noise.freqvec = 1;

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1];
%         noise.scalevec = [5; 1];
%         noise.freqvec = [1; 0.1];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 0; 1];
%         noise.scalevec = [5; 1; 1];
%         noise.freqvec = [1; 0.5; 0.1];
        
%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 0; 1];
%         noise.scalevec = [10; 5; 5];
%         noise.freqvec = [1; 0.5; 0.1];        

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1; 0 ];
%         noise.scalevec = [20; 5; 5 ];
%         noise.freqvec = [1; 5; 0.1 ];

%         % Multi-sinusoid signal
%         noise.cos1_sin0 = [1; 1; 1; 0];
%         noise.scalevec = [10; 1; 5; 5];
%         noise.freqvec = [1; 5; 0.1; 0.5];
        
        % Learning interval [0, t_f]
        tf = 3;
        
        % Time to simulate after learning [t_f, t_f + tsim]
        tsim = 10;
        
        % Number of iterations to terminate after (IRL and RADP only)
        istar = 20;
        
        % Initial stabilizing policy (if algorithm uses it)
%         u_0.tag = 'comp_vrabie_lewis_2009_hard';   
        u_0.tag = 'vrabie_lewis_2009_hard_min_critic'; 
%         u_0.tag = 'vrabie_lewis_2009_hard_x_2_only'; 
%         u_0.tag = 'vrabie_lewis_2009_hard_testing'; 

        % ***********************
        %
        % CRITIC BASIS   
        %  
        
        % Default basis for critic NN; i.e., \hat{V}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        
        basis_critic.tag = 'vrabie_lewis_2009_hard_critic_7terms_no_x_1^2';
%         basis_critic.tag = 'vrabie_lewis_2009_hard_critic_7terms_no_x_2^4';
%         basis_critic.tag = 'vrabie_lewis_2009_hard_critic_8terms';
%         basis_critic.tag = 'vrabie_lewis_2009_hard_critic_nonmin';         
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_debug';
%         basis_critic.tag = 'order_2_degree_4';

%         % All even monomials of total degree <= 4
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'even';   
%         basis_critic.K = 4;
%         basis.critic.IC_type = 'none';

%         % All monomials of total degree = 2
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'total_deg_K';   
%         basis_critic.K = 2;
%         basis.critic.IC_type = 'none';

        % ICs
        basis_critic.c0 = [3/2 ; 3; zeros(5,1)];
        
        % ***********************
        %
        % ACTOR BASIS -- g(x) UNKNOWN
        %  
        % Basis for actor NN without g(x); i.e., \hat{\mu}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        %
        
        basis_actor_no_g.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';   

%         % All even monomials of total degree <= 2
%         basis_actor_no_g.tag = 'monomial';
%         basis_actor_no_g.type = 'even';   
%         basis_actor_no_g.K = 2;
%         basis_actor_no_g.IC_type = 'none';


        % ***********************
        %
        % HAMILTONIAN NETWORK BASIS (VI ONLY)
        %        
        
        % Hamiltonian NN basis \Theta(x)
%         Theta.tag = 'comp_vrabie_lewis_2009_hard_vi_theta';
        Theta.tag = 'comp_vrabie_lewis_2009_hard_vi_theta_w_Q';
        
        % Store parameters
        basis_hamiltonian.Theta = Theta;

        % ***********************
        %
        % OPTIMAL BASES, WEIGHTS
        %        
        
        % Tags corresponding to optimal bases, associated optimal weights
        basis_critic_opt.tag = 'comp_vrabie_lewis_2009_hard_critic_min';
        basis_critic_opt.c_star = [0.5; 1; 1];
        
        basis_actor_no_g_opt.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';
        basis_actor_no_g_opt.w_star = [-1; -2];
        basis_actor_no_g_opt.w_star_VI = -2 * [-1; -2];
        
        basis_hamiltonian_opt.tag = 'comp_vrabie_lewis_2009_hard_vi_theta';
        basis_hamiltonian_opt.v_star = [1; 4; 4;];
        
        % Optimal weights associated with the bases used (NOT with the
        % optimal bases)
%         c_star = [0.5; 1; 1; 0];
        w_star = [-1; -2]; 
        w_star_VI = -2 * w_star; 
        v_star = [1; 4; 4; -1];  
        
        % Store if each of the bases used are indeed the optimal bases
        % (or contain the optimal basis functions as a subset)
        isopt_critic = 0;
        isopt_actor_no_g = 1;
        isopt_hamiltonian = 1;
        
        
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        
        
        % Do plots for each individual preset in the group
        do_individual_plots = 1;
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 1;
                      
        % Whether or not to plot the optimal weight values on surface plots
        % for critic weights
        plot_opt_weight = 1;

        
        % ***********************
        %
        % STORE SETTINGS
        %
        
        group_settings.Q = Q;
        group_settings.R = R;
        group_settings.noise = noise;
        group_settings.tf = tf;
        group_settings.tsim = tsim;
        
        group_settings.istar = istar;

        group_settings.alg_list = alg_list;
        group_settings.numalgs = numalgs;
        group_settings.x0_cell = x0_cell;
        group_settings.numICs = numICs;
        group_settings.x1vec = x1vec;
        group_settings.nx1 = nx1;
        group_settings.x2vec = x2vec;
        group_settings.nx2 = nx2;
        group_settings.xinds = xinds;
        group_settings.remove_0_IC = remove_0_IC;
        group_settings.ind_0_IC = ind_0_IC;
        
        group_settings.plot_opt_weight = plot_opt_weight;
%         group_settings.c_star = c_star;
        group_settings.w_star = w_star;
        group_settings.w_star_VI = w_star_VI;
        group_settings.v_star = v_star;

        group_settings.basis_critic_opt = basis_critic_opt;
        group_settings.basis_actor_no_g_opt = basis_actor_no_g_opt;
        group_settings.basis_hamiltonian_opt = basis_hamiltonian_opt;
        
        group_settings.isopt_critic = isopt_critic;
        group_settings.isopt_actor_no_g = isopt_actor_no_g;
        group_settings.isopt_hamiltonian = isopt_hamiltonian;
        
        group_settings.aux_plot_x0s = aux_plot_x0s;             
                
        
    % *********************************************************************
    %
    % VRABIE, LEWIS (2009) -- NONLINEAR "HARD" EXAMPLE
    %
    
    case 'vrabie_lewis_2009_hard_ex'
        
        %%   
        % ***********************
        %
        % PRESET LIST
        %      
        
        preset_list = {'vrabie_lewis_2009_hard_ex'};

        % ***********************
        %
        % PRESET GROUP SHARED SETTINGS
        %
        % Note: This is a single-preset group. As such, most of the
        % preset settings are declared in config.m        
        %        
        
        % System tag
        sys.tag = 'vrabie_lewis_2009_hard';

        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        

        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 0;        
        
        % Do plots for each individual preset in the group
        do_individual_plots = 1;
        
    % *********************************************************************
    %
    % VAMVOUDAKIS, LEWIS (2010) -- LINEAR F16 EXAMPLE
    %
    
    case 'vamvoudakis_lewis_2010_F16_lin_ex'
        
        %%   
        % ***********************
        %
        % PRESET LIST
        %        
        
        preset_list = {'vamvoudakis_lewis_2010_F16_lin_ex'};

        % ***********************
        %
        % PRESET GROUP SHARED SETTINGS
        %
        % Note: This is a single-preset group. As such, most of the
        % preset settings are declared in config.m        
        %        
        
        % System tag
        sys.tag = 'vamvoudakis_lewis_2010_F16_lin';

        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;    
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 0;
        
        % Do plots for each individual preset in the group
        do_individual_plots = 1;
        
    % *********************************************************************
    %
    % VAMVOUDAKIS, LEWIS (2010) -- NONLINEAR EXAMPLE
    %
    
    case 'vamvoudakis_lewis_2010_nonlin_ex'  
   
        %%   
        % ***********************
        %
        % PRESET LIST
        %        
        
        preset_list = {'vamvoudakis_lewis_2010_nonlin_ex'};

        % ***********************
        %
        % PRESET GROUP SHARED SETTINGS
        %
        % Note: This is a single-preset group. As such, most of the
        % preset settings are declared in config.m        
        %        
        
        % System tag
        sys.tag = 'vamvoudakis_lewis_2010_nonlin';
          
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;

        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 0;        
        
        % Do plots for each individual preset in the group
        do_individual_plots = 1;
        

    % *********************************************************************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE EXAMPLE
    %        
        
    case 'jiang_jiang_2014_engine_ex'
        
        %%   
        % ***********************
        %
        % PRESET LIST
        %        
        
        preset_list = {'jiang_jiang_2014_engine_ex'};

        % ***********************
        %
        % PRESET GROUP SHARED SETTINGS
        %
        % Note: This is a single-preset group. As such, most of the
        % preset settings are declared in config.m        
        %        
        
        % System tag
        sys.tag = 'jiang_jiang_2014_engine';
        
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;     
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 0;
        
        % Do plots for each individual preset in the group
        do_individual_plots = 1;
        
    % *********************************************************************
    %
    % T. BIAN, Z.-P. JIANG (2021) -- NONLINEAR EXAMPLE
    %
    
    case 'bian_jiang_2021_nonlin'
        
        %%   
        % ***********************
        %
        % PRESET LIST
        %        
        
        preset_list = {'bian_jiang_2021_nonlin'};

        % ***********************
        %
        % PRESET GROUP SHARED SETTINGS
        %
        % Note: This is a single-preset group. As such, most of the
        % preset settings are declared in config.m        
        %        
        
        % System tag
        sys.tag = 'bian_jiang_2021_nonlin_ex'; 
        
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 0;
        
        % Do plots for each individual preset in the group
        do_individual_plots = 1;

        
    % *********************************************************************
    %
    % RADP ON 2ND ORDER SYSTEM -- SWEEP x_0
    %
    
    case 'radp_2ndorder_sweep_x0'
        
        %%   
        
        % *****************************************************************
        %
        % IC VECTORS TO SWEEP
        %   
        
        % Vector of x_1(0)
%         x1vec = (-1:0.25:1)';
%         x1vec = (-1:0.5:1)';
        x1vec = [-1;1];
%         x1vec = [-1;0;1];
%         x1vec = [0];
%         x1vec = [-1;0.5;1];        
        
        % Vector of x_2(0)
%         x2vec = (-1:0.25:1)';
%         x2vec = (-1:0.5:1)';
        x2vec = [-1;1];
%         x2vec = [-1;0;1];
%         x2vec = [0];
%         x2vec = [-1;0.75;1];   

        % Size of IC vectors
        nx10 = size(x1vec,1);
        nx20 = size(x2vec,1);
        
        % Replace any value of x_2(0) that is zero with a small number
        inds = x2vec == 0;
        x2vec(inds) = 0.05 * max([max(abs(x2vec)) 1]);
        
        % Store these IC vectors in a cell array for indexing
        x0_cell = cell(nx10,nx20);
        for i = 1:nx10
            for j = 1:nx20
                x0_cell{i,j} = [x1vec(i) ; x2vec(j)];
            end
        end
        
        % Matrix of ICs to plot estimate optimal value, policy for.
        % E.g.,  if plot_opt_x0s(2,:) = [1 -1], then the second set of ICs
        % for which the estimate optimal value, policy will be plotted for
        % will be for the IC x_0 = [1; -1]
        aux_plot_x0s = [
                        -1  -1 
                        -1  1
                        1   -1
                        1   1
                                ];
                            
        % *****************************************************************
        %
        % PRESET LIST
        %      
        
        preset_list = cell(nx10*nx20, 1);
        
        for i = 1:nx10*nx20
            preset_list{i} = 'radp_2ndorder_sweep_x0';
        end

        % *****************************************************************
        %
        % PRESET GROUP SHARED SETTINGS      
        %        

        % System tag
        sys.tag = 'vrabie_lewis_2009_hard';
        
        % State penalty function
        Q = 'vrabie_lewis_2009_hard_ex';
        
        % Control penalty matrix
        R = 1;
        
        % Probing noise signal
        noise.tag = 'radp_2ndorder_sweep_x0';
        
        % Learning interval [0, t_f]
        tf = 10;
        
        % Time to simulate after learning [t_f, t_f + tsim]
        tsim = 10;
        
        % Number of iterations to perform before terminating
        istar = 25;

        % ***********************
        %
        % CRITIC BASIS   
        %  
        
        % Default basis for critic NN; i.e., \hat{V}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        
        basis_critic.tag = 'order_2_degree_4';
        
%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_min';

%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_debug';

%         % All even monomials of total degree <= 4
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'even';   
%         basis_critic.K = 4;
%         basis_critic.IC_type = 'none';

%         % All monomials of total degree = 2
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'total_deg_K';   
%         basis_critic.K = 2;
%         basis_critic.IC_type = 'none';
        
        % ***********************
        %
        % ACTOR BASIS -- g(x) UNKNOWN
        %  
        
        % Basis for actor NN without g(x); i.e., \hat{\mu}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        
        basis_actor_no_g.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';   

%         % All even monomials of total degree <= 2
%         basis_actor_no_g.tag = 'monomial';
%         basis_actor_no_g.type = 'even';   
%         basis_actor_no_g.K = 2;
%         basis_actor_no_g.IC_type = 'none';
   
        
        % Initial stabilizing policy (if algorithm uses it)
        u_0.tag = 'comp_vrabie_lewis_2009_hard';
                      
        % Whether or not to plot the optimal weight values on surface plots
        % for critic weights
        plot_opt_weight = 1;
        
        % Optimal weights
        c_star = [0.5; 0; 1; 0; 0; 0; 0; 1];
        w_star = [-1; -2];        

        
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 1;
        
        % Do plots for each individual preset in the group
        do_individual_plots = 0;
        
        % ***********************
        %
        % STORE SETTINGS
        %
        
        group_settings.Q = Q;
        group_settings.R = R;
        group_settings.noise = noise;
        group_settings.tf = tf;
        group_settings.tsim = tsim;
%         group_settings.x0 = x0;
        group_settings.istar = istar;

        group_settings.x0_cell = x0_cell;
        group_settings.x1vec = x1vec;
        group_settings.nx10 = nx10;
        group_settings.x2vec = x2vec;
        group_settings.nx20 = nx20;
        
        group_settings.plot_opt_weight = plot_opt_weight;
        group_settings.c_star = c_star;
        group_settings.w_star = w_star;
        
        group_settings.plot_opt_x0s = aux_plot_x0s;

        
    % *********************************************************************
    %
    % RADP ON CART INVERTED PENDULUM SYSTEM -- SWEEP x_0
    %
    
    case 'radp_cip_sweep_x0'
        
        %%   
        
        % *****************************************************************
        %
        % IC VECTORS TO SWEEP
        %   
        
        % Vector of x(0)
%         xvec = (-1:0.25:1)';
%         xvec = (-1:0.5:1)';
%         xvec = [-1;1];
         xvec = [-0.5;0.5];
%         xvec = [-1;0;1];
%         xvec = [0];
%         xvec = [-1;0.5;1];        
        
        % Vector of x_2(0)
%         thetavec = (-1:0.25:1)';
%         thetavec = (-1:0.5:1)';
        thetavec = [-0.5;0.5];
%         thetavec = [-1;0;1];
%         thetavec = [0];
%         thetavec = [-1;0.75;1];   

        % Size of IC vectors
        nx0 = size(xvec,1);
        ntheta0 = size(thetavec,1);
        
%         % Replace any value of x_2(0) that is zero with a small number
%         inds = xthetavec == 0;
%         xthetavec(inds) = 0.05 * max([max(abs(xthetavec)) 1]);
        
        % Store these IC vectors in a cell array for indexing
        x0_cell = cell(nx0,ntheta0);
        for i = 1:nx0
            for j = 1:ntheta0
                x0_cell{i,j} = [xvec(i) ; thetavec(j)];
            end
        end
        
        % Matrix of ICs to plot estimate optimal value, policy for.
        % E.g.,  if plot_opt_x0s(2,:) = [1 -1], then the second set of ICs
        % for which the estimate optimal value, policy will be plotted for
        % will be for the IC x_0 = [1; -1]
        aux_plot_x0s = [
                        -0.5  -0.5 
                        -0.5  0.5
                        0.5   -0.5
                        0.5   0.5
                                ];
        
        % *****************************************************************
        %
        % PRESET LIST
        %      
        
        preset_list = cell(nx0*ntheta0, 1);
        
        for i = 1:nx0*ntheta0
            preset_list{i} = 'radp_cip_sweep_x0';
        end

        % *****************************************************************
        %
        % PRESET GROUP SHARED SETTINGS      
        %        

        % System tag
        sys.tag = 'cip';
        
        % State penalty function
        Q = eye(4);
        
        % Control penalty matrix
        R = 1;
        
        % Probing noise signal
        noise.tag = 'radp_cip_sweep_x0';
        
        % Learning interval [0, t_f]
        tf = 10;
        
        % Time to simulate after learning [t_f, t_f + tsim]
        tsim = 10;

        % ***********************
        %
        % CRITIC BASIS   
        %  
        
        % Default basis for critic NN; i.e., \hat{V}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        
%         % All even monomials of total degree <= 4
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'even';   
%         basis_critic.K = 4;
%         basis_critic.IC_type = 'none'; 

        % All monomials of total degree = 2
        basis_critic.tag = 'monomial';
        basis_critic.type = 'total_deg_K';   
        basis_critic.K = 2;
        basis_critic.IC_type = 'none'; 
        
        % ***********************
        %
        % ACTOR BASIS -- g(x) UNKNOWN
        %  
        
        % Basis for actor NN without g(x); i.e., \hat{\mu}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        

        % All odd monomials of total degree <= 3
        basis_actor_no_g.tag = 'monomial';
        basis_actor_no_g.type = 'odd';   
        basis_actor_no_g.K = 3;
        basis_actor_no_g.IC_type = 'none'; 

%         % Manually-chosen odd monomials of total degree <= 3
%         basis_actor_no_g.tag = 'monomial';
%         basis_actor_no_g.type = 'custom';  
%         basis_actor_no_g.IC_type = 'custom';
%         basis_actor_no_g.dmat = [   
%                                     1 0 0 0     % x_1
%                                     0 1 0 0     % x_2
%                                     0 0 1 0     % x_3
%                                     0 0 0 1     % x_4
% %                                     0 0 3 0     % x_3^3
% %                                     2 0 1 0     % x_1^2 * x_3
% %                                     1 0 2 0     % x_1 * x_3^2
% %                                     0 0 2 1     % x_3^2 * x_4
% %                                     1 1 1 0     % x_1 * x_2 * x_3
% %                                     1 0 1 1     % x_1 * x_3 * x_4
%                                             ];        
        
        
        % Initial stabilizing policy (if algorithm uses it)
        u_0.tag = 'comp_cip';

        
        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 1;
        
        % Do plots for each individual preset in the group
        do_individual_plots = 0;    
        

        % ***********************
        %
        % STORE SETTINGS
        %          
        
        group_settings.Q = Q;
        group_settings.R = R;
        group_settings.noise = noise;
        group_settings.tf = tf;
        group_settings.tsim = tsim;
%         group_settings.x0 = x0;

        group_settings.x0_cell = x0_cell;
        group_settings.xvec = xvec;
        group_settings.nx0 = nx0;
        group_settings.thetavec = thetavec;
        group_settings.ntheta0 = ntheta0;
        
        group_settings.plot_opt_x0s = aux_plot_x0s;
    
        
    % *********************************************************************
    %
    % ALL METHOD COMPARISON OF D. VRABIE, F.L. LEWIS (2009) NONLINEAR
    % "HARD" EXAMPLE (SEC. 6.2)
    %
    
    case 'comp_vrabie_lewis_2009_hard'
        
        %%   
        % ***********************
        %
        % PRESET LIST
        %      
        
        % DEBUGGING: One method at a time
%         preset_list = {'comp_vrabie_lewis_2009_hard_irl'};
%         preset_list = {'comp_vrabie_lewis_2009_hard_spi'};
%         preset_list = {'comp_vrabie_lewis_2009_hard_radp_matched'};
        preset_list = {'comp_vrabie_lewis_2009_hard_vi'};
        
%         preset_list = {'comp_vrabie_lewis_2009_hard_irl';...
%                     'comp_vrabie_lewis_2009_hard_spi';...
%                     'comp_vrabie_lewis_2009_hard_radp_matched';...
%                     'comp_vrabie_lewis_2009_hard_vi';...
%                         };

        % ***********************
        %
        % PRESET GROUP SHARED SETTINGS      
        %        
        
        % System tag
        sys.tag = 'vrabie_lewis_2009_hard';
        
        % State penalty function
        Q = 'vrabie_lewis_2009_hard_ex';
        
        % Control penalty matrix
        R = 1;
        
        % Probing noise signal
%         noise.tag = 'comp_vrabie_lewis_2009_hard';

        % See eval_noise.m for documentation of how the 'sum_sinusoids'
        % exploration noise option works.
        noise.tag = 'sum_sinusoids';
%         % A cos(t)
%         noise.cos1_sin0 = 1;
%         noise.scalevec = 5;
%         noise.freqvec = 1;       
        % A sin(t)
        noise.cos1_sin0 = 0;
        noise.scalevec = 5;
        noise.freqvec = 1;  
        
        % Learning interval [0, t_f]
        tf = 10;
        
        % Time to simulate after learning [t_f, t_f + tsim]
        tsim = 10;

        % ***********************
        %
        % CRITIC BASIS   
        %  
        
        % Default basis for critic NN; i.e., \hat{V}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
%         basis_critic.tag = 'order_2_degree_4';
        
        basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_min';

%         basis_critic.tag = 'comp_vrabie_lewis_2009_hard_critic_debug';

%         % All even monomials of total degree <= 4
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'even';   
%         basis_critic.K = 4;
%         basis_critic.IC_type = 'none';

%         % All monomials of total degree = 2
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'total_deg_K';   
%         basis_critic.K = 2;
%         basis_critic.IC_type = 'none';
        
        % ***********************
        %
        % ACTOR BASIS -- g(x) UNKNOWN
        %  
        
        % Basis for actor NN without g(x); i.e., \hat{\mu}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
        basis_actor_no_g.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';   


%         % All even monomials of total degree <= K
%         basis_actor_no_g.tag = 'monomial';
%         basis_actor_no_g.type = 'even';   
%         basis_actor_no_g.K = 4;
%         basis_actor_no_g.IC_type = 'none';
        
        % Initial stabilizing policy (if algorithm uses it)
%         u_0.tag = 'comp_vrabie_lewis_2009_hard';
        u_0.tag = 'vrabie_lewis_2009_hard_min_critic';
               
        % Initial conditions
        x0 = [1; 1];
%         x0 = [-2; 1];

        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;        
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 0;
        
        % Do plots for each individual preset in the group
        do_individual_plots = 1;
        
        
        % ***********************
        %
        % STORE SETTINGS
        %
        
        group_settings.Q = Q;
        group_settings.R = R;
        group_settings.noise = noise;
        group_settings.tf = tf;
        group_settings.tsim = tsim;
        group_settings.x0 = x0;

        
    % *********************************************************************
    %
    % ALL METHOD COMPARISON OF CART INVERTED PENDULUM SYSTEM
    %
    
    case 'comp_cip'
        
        %%   
        % *****************************************************************
        %
        % PRESET LIST
        %      
        
        % DEBUGGING: One method at a time
        preset_list = {'comp_cip_irl'};
%         preset_list = {'comp_cip_spi'};
%         preset_list = {'comp_cip_radp_matched'};
%         preset_list = {'comp_cip_vi'};
%         preset_list = {'comp_cip_test'};
        
%         preset_list = {'comp_cip_irl';...
%                     'comp_cip_spi';,...
%                     'comp_cip_radp_matched';...
%                     'comp_cip_vi';...
%                         };

        % *****************************************************************
        %
        % PRESET GROUP SHARED SETTINGS      
        %        
        
        % System tag
        sys.tag = 'cip';
        
        % State penalty function
        Q = eye(4);
        
        % Control penalty matrix
%         R = 0.1;
        R = 1;
        
        % Probing noise signal
        noise.tag = 'comp_cip';
        
        % Learning interval [0, t_f]
%         tf = 1;
%         tf = 5;
        tf = 10;
%         tf = 20;
%         tf = 100;
        
        % Time to simulate after learning [t_f, t_f + tsim]
%         tsim = 0.1;
        tsim = 10;

        % Obtained from LQR on linearization
        u_0_K = [ -31.6228;  -50.2741; -247.4658;  -63.6293];        
        
        % ***********************
        %
        % CRITIC BASIS   
        %  
        
        % Default basis for critic NN; i.e., \hat{V}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
%         basis_critic.tag = 'order_2_degree_4';

%         % Even monomials of total degree <= 4
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'even';   
%         basis_critic.K = 4;
%         basis_critic.IC_type = 'none';
        
%         % All monomials of total degree <= 2
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'all';   
%         basis_critic.K = 2;     
%         basis_critic.IC_type = 'none';
        
        % Monomials of total degree = 2
        basis_critic.tag = 'monomial';
        basis_critic.type = 'total_deg_K';   
        basis_critic.K = 2;
        basis_critic.IC_type = 'none';
        

%         % Manually-chosen even monomials of total degree =2
%         basis_critic.tag = 'monomial';
%         basis_critic.type = 'custom';   
%         basis_critic.IC_type = 'custom';
%         basis_critic.dmat = [   
%                                     2 0 0 0     % x_1^2
%                                     0 2 0 0     % x_2^2
%                                     0 0 2 0     % x_3^2
%                                     0 0 0 2     % x_4
% %                                     1 1 0 0     % x_1 * x_2
%                                     1 0 1 0     % x_1 * x_3
% %                                     1 0 0 1     % x_1 * x_4
%                                     0 1 1 0     % x_2 * x_3
% %                                     0 1 0 1     % x_2 * x_4
% %                                     0 0 1 1     % x_3 * x_4
%                                             ];
  

%         % Initial conditions for actor with known g(x)
%         % The rows of 'dmat_ic' correspond to the monomials in the critic
%         % 
%         % basis activation functions \Phi(x) which, once partially
%         % differentiated for the actor, yield the desired rows of \nabla
%         % \Phi(x) to set with the ICs in c0_ic. E.g., dmat_ic(1,:)
%         % corresponds to x_1^2, which when partially differentiated yields
%         % [2*x_1 0 0 0] in \nabla \Phi(x). c0_ic(1) contains the desired
%         % actor IC to weight this row of \nabla \Phi(x) with.
%         %
%         dmat_ic = [ 
%                     2 0 0 0 
%                     0 2 0 0
%                     0 0 2 0
%                     0 0 0 2 
%                             ];
%         c0_ic = 1 * [1 ; 10 ; 1 ; 10];
%         
%         % Store these parameters
%         basis_actor_g_known.dmat_ic = dmat_ic;
%         basis_actor_g_known.c0_ic = c0_ic;

        % This gain matrix K stabilizes the system:
        %
        % (A, -1/2 B R^{-1} B^T)
        %
        % Where (A, B) are the linearization of the nonlinear system.
        K_mod = [   -0.0000   -0.0000   -0.0000   -0.0000
                    0.4472    0.9692   12.4134    2.7487
                    0         0         0         0
                    -0.8944   -1.9385  -24.8269   -5.4973   ];
        
        n = 4;
        dmat_ic = zeros(n^2,4);   
        c0_ic = zeros(n^2,1);
        count = 1;
        for i = 1:n
            for j = 1:n
                % Extaract the gain for this entry
                k_ij = K_mod(i,j);
                % Look for the activation function x_i * x_j
                degs = zeros(1,n);
                degs(i) = 1;
                degs(j) = 1;
                dmat_ic(count,:) = degs;
                c0_ic(count) = k_ij;
                count = count + 1;
            end
        end
            
        % Store these parameters
        basis_critic.dmat_ic = dmat_ic;
        basis_critic.c0_ic = c0_ic;      

        
        % ***********************
        %
        % ACTOR BASIS -- g(x) UNKNOWN
        %  
        
        % Basis for actor NN without g(x); i.e., \hat{\mu}(x) =
        % \sum_{i=1}^{N} w_i * \phi_i(x).
%         basis_actor_no_g.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';

        % Odd monomials of total degree <= 3
        basis_actor_no_g.tag = 'monomial';
        basis_actor_no_g.type = 'odd';   
        basis_actor_no_g.K = 3;
        basis_actor_no_g.IC_type = 'none';
        
%         % All monomials of total degree <= 3
%         basis_actor_no_g.tag = 'monomial';
%         basis_actor_no_g.type = 'all';   
%         basis_actor_no_g.K = 3;   
%         basis_actor_no_g.IC_type = 'none';
        
%         % Manually-chosen odd monomials of total degree <= 3
%         basis_actor_no_g.tag = 'monomial';
%         basis_actor_no_g.type = 'custom';
%         basis_actor_no_g.IC_type = 'custom';
%         basis_actor_no_g.dmat = [   
%                                     1 0 0 0     % x_1
%                                     0 1 0 0     % x_2
%                                     0 0 1 0     % x_3
%                                     0 0 0 1     % x_4
% %                                     0 0 3 0     % x_3^3
% %                                     2 0 1 0     % x_1^2 * x_3
% %                                     1 0 2 0     % x_1 * x_3^2
% %                                     0 0 2 1     % x_3^2 * x_4
% %                                     1 1 1 0     % x_1 * x_2 * x_3
% %                                     1 0 1 1     % x_1 * x_3 * x_4
%                                             ];
%         % ICs for manually-chosen basis
%         dmat_ic = [ 1 0 0 0 
%                     0 1 0 0
%                     0 0 1 0
%                     0 0 0 1 ];
%         c0_ic = - u_0_K;
%         % Store these parameters
%         basis_actor_no_g.dmat_ic = dmat_ic;
%         basis_actor_no_g.c0_ic = c0_ic;
        
      
        % Initial stabilizing policy (if algorithm uses it)
        u_0.tag = 'comp_cip';
        
               
        % Initial conditions
%         x0 = [1; 0.1; 0.5; 0];
%         x0 = [1; -0.1; -0.5; 0];
        x0 = [-0.5; 0; -0.5; 0];

        % ***********************
        %
        % PLOT SETTINGS
        %
        
        % Include legend in plots
        dolegend = 0;   
        
        % Bool which indicates if this preset group is an IC sweep
        is_x0_sweep = 0;
        
        % Do plots for each individual preset in the group
        do_individual_plots = 1;        
        
        % ***********************
        %
        % STORE SETTINGS
        %
        
        group_settings.Q = Q;
        group_settings.R = R;
        group_settings.noise = noise;
        group_settings.tf = tf;
        group_settings.tsim = tsim;
        group_settings.x0 = x0;
                        
        
        
    % *********************************************************************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %   
    
    otherwise
        
        error('*** ERROR: PRESET GROUP TAG NOT RECOGNIZED ***');  
       
end    


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INITIALIZE SYSTEM PARAMETERS AND PLOT SETTINGS
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% Configure settings
[sys, sys_plot_settings] = config_sys(sys);

% Store system parameters
group_settings.sys = sys;

% Store plot settings
group_settings.sys_plot_settings = sys_plot_settings;


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INITIALIZE BASES
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************


% List of networks by variable name
NN_list = {    
                'basis_critic'
                'basis_actor_no_g'
%                 'basis_actor_g_known'
                'basis_hamiltonian' 
                                        };

% Number of network types
num_NN = size(NN_list, 1);
                                    
% Basis initialization settings
b_sett.sys = sys;
b_sett.alg = 'debug';
b_sett.Q = Q;
b_sett.R = R;
                                    
% ***********************
%
% INITIALIZE BASES
%

for i = 1:num_NN
   
    % Extract current network type
    NN_i = NN_list{i};
    
    % Check if the current network type exists as a variable
    if exist(NN_i, 'var')
        
        % Extract the basis
        switch NN_i
            case 'basis_critic'
                basis_i = basis_critic;
            case 'basis_actor_no_g'
                basis_i = basis_actor_no_g;
            case 'basis_actor_g_known'
                basis_i = basis_critic;
            case 'basis_hamiltonian'
                basis_i = basis_hamiltonian.Theta;
        end
        
        % Initialize the basis parameters
        basis_i = config_basis(basis_i, b_sett);
        
        % Store the initialized basis
        switch NN_i
            case 'basis_critic'
                basis_critic = basis_i;
                group_settings.basis_critic = basis_critic;
            case 'basis_actor_no_g'
                basis_actor_no_g = basis_i;
                group_settings.basis_actor_no_g = basis_actor_no_g;
            case 'basis_actor_g_known'
                basis_critic = basis_i;
                group_settings.basis_actor_g_known = basis_critic;
            case 'basis_hamiltonian'
                basis_hamiltonian.Theta = basis_i;
                group_settings.basis_hamiltonian.Theta = basis_hamiltonian.Theta;
        end
        
    end
    
end


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INITIALIZE INITIAL STABILIZING POLICY PARAMETERS (IF REQUIRED)
%
% See eval_u.m for details.
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************


if exist('u_0', 'var')
   
    % Check if additional fields need to be declared for the initial
    % stabilizing policy type 
    switch u_0.tag
        case 'actor_g_known'
            u_0.basis = basis_critic;
            u_0.sys = sys;            
            u_0.R = R;
        case 'actor_g_known_lqr'
            u_0.basis = basis_critic;
            u_0.sys = sys;
            u_0.R = R;
        otherwise
            % Do nothing
    end
    
    % Store the initial policy
    group_settings.u_0 = u_0;
    
end



%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% STORE GROUP SETTINGS WHICH ARE ALWAYS DECLARED
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% Preset group
group_settings.preset_group = preset_group;

% Do plots for each individual preset in the group
group_settings.do_individual_plots = do_individual_plots;

% Bool which indicates if this preset group is an IC sweep
group_settings.is_x0_sweep = is_x0_sweep;

% Include legend in plots
group_settings.dolegend = dolegend; 
