% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% COMPARISON OF CT RL-BASED OPTIMAL CONTROL TECHNIQUES
%
% Brent Wallace  
%
% 2021-11-06
%
% Main file for comparison of CT RL-based optimal control techniques.
% Considers infinite-horizon, input-affine nonlinear systems.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% 
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CONFIG
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% ***********************
%
% CLEAR VARIABLES, COMMAND WINDOW, FIGURES
%
clear
clc
close all

% ***********************
%
% INCLUDE UTILITY FOLDER
%
addpath('util');

% ***********************
%
% INCLUDE CONFIG FUNCTIONS FOLDER
%
addpath('config');

% ***********************
%
% INCLUDE EVALUATION FUNCTIONS FOLDER
%
addpath('eval_functs');

% ***********************
%
% INCLUDE ALGORITHM FOLDER
%
addpath('algs');

% ***********************
%
% INCLUDE PLOT FUNCTIONS FOLDER
%
addpath('plot');

% ***********************
%
% INCLUDE DATA MANAGEMENT FUNCTIONS FOLDER
%
addpath('data_mgmt');

% ***********************
%
% FIGURES
%
savefigs = 1;               % Save figures control
relpath = '00 figures/';  	% Relative file path for saved figures

% *************************************************************************
%
% METHOD/SYSTEM/DESIGN PRESET GROUPS
%
% These tags correspond to the group of presets to be executed. Each preset
% within the group contains the specific
%
%   Algorithm/methodology (e.g., IRL, SPI, etc.)
%   System (e.g., Vrabie, Lewis (2009) 2nd order example system)
%   Design (with numerical design parameters)
%
% NOTE: See config.m for a complete list of algorithm, system, and preset
% options.
%
% Preset group options:
%
%   CS_2ndorder_N1_3
%       Case study I. All methodologies implemented on the 2nd order
%       academic example in D. Vrabie, F.L. Lewis (2009), Sec. 6.2. Minimal
%       critic, actor, and Hamiltonian bases. Sweep ICs x0.
%
%   CS_2ndorder_N1_4
%       Case study II. All methodologies implemented on the 2nd order
%       academic example in D. Vrabie, F.L. Lewis (2009), Sec. 6.2. Minimal
%       bases, except critic with N_1 = 4 terms (activation function x_1 *
%       x_2 added). Sweep ICs x0.
%
%   CS_2ndorder_N1_7_no_x2p4
%       Case study III. All methodologies implemented on the 2nd order
%       academic example in D. Vrabie, F.L. Lewis (2009), Sec. 6.2. Minimal
%       bases, except critic with N_1 = 7 terms. Essential activation
%       function x_2^4 removed. Sweep ICs x0.
%
%   CS_2ndorder_N1_8
%       Case study. All methodologies implemented on the 2nd order
%       academic example in D. Vrabie, F.L. Lewis (2009), Sec. 6.2. Minimal
%       bases, except critic with N_1 = 8 terms. Sweep ICs x0.
%
%   CS_2ndorder_N1_7_no_x1p2
%       Case study. All methodologies implemented on the 2nd order
%       academic example in D. Vrabie, F.L. Lewis (2009), Sec. 6.2. Minimal
%       bases, except critic with N_1 = 7 terms. Essential activation
%       function x_1^2 removed. Sweep ICs x0.
%
%   CS_cip
%       Case study. All methodologies implemented on the cart inverted
%       pendulum system. Sweep ICs x0.
%
%   vrabie_lewis_2009_hard_ex
%       (Single-design group) Implementation of the example in Sec. 6.2
%       (see preset description for further detail).
%
%   vamvoudakis_lewis_2010_F16_lin_ex
%       (Single-design group) Implementation of the example in Sec. 5.2
%       (see preset description for further detail).
%
%   vamvoudakis_lewis_2010_nonlin_ex
%       (Single-design group) Implementation of the example in Sec. 5.2
%       (see preset description for further detail).
%
%   jiang_jiang_2014_engine_ex
%       (Single-design group) Implementation of the example in Sec. V. A.
%       (see preset description for further detail).
%
%   bian_jiang_2021_nonlin
%       (Single-design group) Implementation of the example in Sec. V. A.
%       (see preset description for further detail).
%
%   radp_2ndorder_sweep_x0
%       RADP with matched uncertainty implemented on the 2nd order academic
%       example in D. Vrabie, F.L. Lewis (2009), Sec. 6.2. Sweep ICs x0.
%
%   radp_cip_sweep_x0
%       RADP with matched uncertainty implemented on the cart inverted
%       pendulum system. Sweep ICs x0.
%
%   comp_vrabie_lewis_2009_hard
%       All methodologies implemented on the example in Sec. 6.2 of D.
%       Vrabie, F.L. Lewis (2009).
%
%   comp_cip
%       All methodologies implemented a cart inverted pendulum example.
%
% *************************************************************************

% preset_group = 'CS_2ndorder_N1_3';
% preset_group = 'CS_2ndorder_N1_4';
% preset_group = 'CS_2ndorder_N1_7_no_x2p4';
preset_group = 'CS_cip';

% preset_group = 'CS_2ndorder_N1_8';
% preset_group = 'CS_2ndorder_N1_7_no_x1p2';

% preset_group = 'vrabie_lewis_2009_hard_ex';
% preset_group = 'vamvoudakis_lewis_2010_F16_lin_ex';
% preset_group = 'vamvoudakis_lewis_2010_nonlin_ex';
% preset_group = 'jiang_jiang_2014_engine_ex';
% preset_group = 'bian_jiang_2021_nonlin';
% 
% preset_group = 'radp_2ndorder_sweep_x0';
% preset_group = 'radp_cip_sweep_x0';

% preset_group = 'comp_vrabie_lewis_2009_hard';
% preset_group = 'comp_cip';



% *************************************************************************
% *************************************************************************
%
% METHOD/SYSTEM/DESIGN PRESET LIST AND SETTINGS
%
% This is a list of tags correspond to the presets to execute for the
% selected preset group. Each preset consists of settings which include the
% specific
%
%   Algorithm/methodology (e.g., IRL, SPI, etc.)
%   System (e.g., Vrabie, Lewis (2009) 2nd order example system)
%   Design (with numerical design parameters)
%
% Below is the selection of the preset list, along with initialization of
% any preset group settings which are shared among all the designs in the
% group (e.g., if all designs share the same system, or the same Q, R
% matrices, etc.)
%
% NOTE: See config.m for a complete list of algorithm, system, and preset
% options.
%
% *************************************************************************
% *************************************************************************

[preset_list, group_settings] = config_preset_group(preset_group);



%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% MAIN LOOP
%
% For each preset in the preset group:
%
%   Run preset configuration to initialize design parameters
%   Run respective preset algorithm and collect algorithm output data 
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% ***********************
%
% INITIALIZATION
%  

% Number of presets to execute in the group
numpresets = size(preset_list, 1);

% ***********************
%
% STORAGE
%  

% Holds algorithm settings for each preset in the group
alg_settings_cell = cell(numpresets, 1);

% Holds algorithm output data for each preset in the group
out_data_cell = cell(numpresets, 1);

for i = 1:numpresets
    
    % Display algorithm count
    disp('************************')
    disp('*')
    disp(['* EXECUTING PRESET    ' num2str(i)...
            '   OUT OF      ' num2str(numpresets)])
    disp('*')
    disp('************************')
    
    % Store the current preset count in group_settings, in case it is
    % needed in the config
    group_settings.presetcount = i;
    
    % *********************************************************************
    %
    % SELECT ALGORITHM, SYSTEM, DESIGN PARAMETERS BASED ON PRESET
    % 
    
    alg_settings = config_preset(preset_list{i}, group_settings);
    
    % *********************************************************************
    %
    % RUN ALGORITHM
    % 
    
    % Get algorithm elapsed time -- start timer
    tic;
    
    switch alg_settings.alg

        % ***********************
        %
        % IRL
        %

        case 'irl'

            out_data = alg_irl(alg_settings);

        % ***********************
        %
        % SPI
        %

        case 'spi'

            out_data = alg_spi(alg_settings);       

        % ***********************
        %
        % RADP (MATCHED UNCERTAINTY)
        %

        case 'radp_matched'

            out_data = alg_radp_matched(alg_settings);             
            
        % ***********************
        %
        % RADP (UNMATCHED UNCERTAINTY)
        %

        case 'radp_unmatched'

            out_data = alg_radp_unmatched(alg_settings);        

        % ***********************
        %
        % VI
        %

        case 'vi'

            out_data = alg_vi(alg_settings);          
            
            
        % ***********************
        %
        % THROW ERROR IF TAG DOES NOT COME UP A MATCH
        %   

        otherwise

            error('*** ERROR: ALGORITHM TAG NOT RECOGNIZED ***');             
        
    end
    
    % Get algorithm elapsed time -- stop timer
    out_data.runtime = toc;

    % *********************************************************************
    %
    % STORAGE
    % 
    
    % Store preset settings
    alg_settings_cell{i} = alg_settings;
    
    % Store algorithm output data
    out_data_cell{i} = out_data;
    
    
end



% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% PLOT/SAVE FIGURES
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% Plot settings
group_settings.savefigs = savefigs;
group_settings.relpath = relpath;


% Call plot function
plot_main(alg_settings_cell, out_data_cell, group_settings);



