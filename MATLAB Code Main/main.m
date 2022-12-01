% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% COMPARISON OF CT RL-BASED OPTIMAL CONTROL TECHNIQUES
%
% Brent Wallace  
%
% 2022-07-30
%
% Main file for comparison of CT RL-based optimal control techniques.
% Considers infinite-horizon, input-affine nonlinear systems.
%
% Associated with the work:
%
% B. A. Wallace, J. Si, "Continuous-Time Reinforcement Learning Control: A
% Review of Theoretical Results, Insights on Performance, and Needs for New
% Designs" TNNLS 2022.
%
%
% ***** PROGRAM EXECUTION FLOW:
%
% * Select preset group (see below for "prest group" definition, options)
% * Configure preset group -- config_preset_group.m
%   * Configure system params - config_sys.m
%   * Initialize NN bases -- config_basis.m
% * For each preset in group:
%   * Configure individual preset's alg hyperparams -- config_preset.m
%   * Run preset's alg -- alg_irl.m, alg_spi.m, alg_radp_matched.m,
%   alg_vi.m
%       * Learning over [0,t_f] -- programs in 'eval_functs' folder
% * Plot/save figures -- plot_main.m
%   * Save alg data to directory "00 figures/...../data/"
%   * State trajectory, control signal plots -- plot_x_u.m
%   * Algorithm-specific plots -- plot_irl.m, plot_spi.m,
%   plot_radp_matched.m, plot_vi.m
%   * (x_0 sweeps only) Plots for x_0 sweeps, sweep statistical data --
%   plot_x0_sweep.m
%   * Miscellaneous plots -- plot_misc.m
%   
% ***** DATA MANAGEMENT:
%
% Figures created by this program are written to the folder "00 figures/"
% in a time-stamped subfolder. Raw program data is written to the same
% location. The program data falls under the following three structures:
%
%   group_settings      (Struct) Contains shared settings common to all
%                       presets (e.g., state penalty matrix Q). In the case
%                       of an x_0 sweep, this will contain the ICs chosen.
%                       This struct is mainly initialized in the program
%                       config_preset_group.m
%   alg_settings_cell   ('numpresets' x 1 cell) Cell array containing the
%                       algorithm settings of each of the presets in the
%                       group. This cell array is mainly initialized in the
%                       program config_preset.m. Each entry's parameters
%                       are algorithm-specific. See alg_irl.m, alg_spi.m,
%                       alg_radp_matched.m, alg_vi.m for detailed
%                       descriptions of each algorithm's hyperparameters
%   out_data_cell       ('numpresets' x 1 cell) After each preset is run,
%                       its algorithm output data (e.g., state trajectory
%                       x(t), NN weight responses, etc.) are stored in this
%                       cell array. The data stored is algorithm-specific.
%                       See alg_irl.m, alg_spi.m, alg_radp_matched.m,
%                       alg_vi.m for detailed descriptions of each
%                       algorithm's output data
%
% After data is generated, the data from two different executions of the
% same preset group may be merged into one new group (e.g., if an x_0 sweep
% is run and then you tweak a hyperparameter for SPI, you can merge the new
% SPI data with the old IRL, RADP, CT-VI data without running these
% algorithms again). See the 'data_mgmt' folder for details, in particular
% merge_x0_sweeps.m for x_0 sweep data management.
%
%
% ***** PLOTTING PREVIOUSLY-GENERATED DATA:
%
% The program does not need to re-run all of the algorithms to generate
% their data and generate plots again. Given data saved in pre-generated
% 'group_settings', 'alg_settings_cell', and 'out_data_cell' structs, the
% function re_plot.m may be used to plot the pre-existing data. See
% re_plot.m for details.
%
% ***** PLOT FORMATTING:
%
% The formatting used for the plots is set programattically by the program
% plot_format.m. The procedure followed to generate plots usually goes as
% follows:
%   Generate plot (title, axes labels, data, etc.) 
%   Specify the figure count in the 'figcount' variable.
%   Format the plot via the following call:
%           p_sett.figcount = figcount;
%           plot_format(p_sett, group_settings);
%       Note: Default formatting options can be over-written for an
%       individual plot by adding additional specs to the 'p_sett' struct
%       (see description of plot_format.m for how to do this).
%   
%
% ***** GENERAL TIPS:
%
% * NOTE: In order to execute properly, the programs must be kept in the
% same file path locations relative to main.m.
% * Set 'savefigs' = 1 (below) to save figures and algorithm data. Data
% gets saved to the relative path specified by the 'relpath' variable
% * For a given evaluation study (e.g., 2nd order system, minimal bases) to
% change which exploration noise is used, or which algorithms are executed
% in the x_0 sweep, go to the respective preset group in
% config_preset_group.m (e.g., 'CS_2ndorder_N1_3') and modify the 'noise'
% and 'alg_list' variables, respectively.
% * Generally speaking, default hyperparameter values are written in
% config_preset_group.m (e.g., learning time t_f, etc.). However, if an
% algorithm does not use the default value of a hyperparameter, its value
% will be initialized individually in config_preset.m (e.g., for eval study
% I, the default learning time is t_f = 10s for RADP, CT-VI. But IRL uses
% t_f = 5s, SPI t_f = 500s).
% * Each algorithm has specific plots (e.g., weight responses, etc.) which
% can be plotted by setting the variable 'do_individual_plots' = 1 in
% group_settings.m (i.e., under the respective preset group). Usually this
% feature is not used, since for an x_0 sweep it will plot weight responses
% for each x_0, each algorithm. To see individual plots for an IC of
% interest, mark it in group_settings.m under the 'aux_plot_x0s' variable.
% Note that the x_0 of interest must be in the sweep, otherwise the program
% will fail.
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
% FIGURES
%
savefigs = 0;               % Save figures control
relpath = '00 figures/';  	% Relative file path for saved figures

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
%   (x_0 sweeps only) IC x_0
%
% The preset group construct is meant to be modular. It can be used to
% execute a single preset (e.g., IRL on an example system for one IC) or
% multiple presets (e.g., IRL, SPI on a sweep of ICs, or IRL with sweeping
% learning time t_f over a vector of values). In an x_0 sweep, each value
% of x_0 corresponds to one preset.
%
% NOTE: See config_preset.m for a complete list of algorithm, system
% options.
%
% ***** PRESET GROUP OPTIONS:
%
% PRESETS FROM PAPER
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
%   CS_cip
%       Case study IV. All methodologies implemented on the cart inverted
%       pendulum system. Sweep ICs x0.
%       *** NOTE: All algorithms break down for this case study. Only the
%       IC [1 0 15deg 0] is tested.
%
% CASE STUDY PRESETS NOT USED IN PAPER
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
% PRESETS DUPLICATING CT-RL ALGORITHM WORK EXAMPLES
%
%   vrabie_lewis_2009_hard_ex
%       (Single-design group) Implementation of the example in Sec. 6.2
%
%   vamvoudakis_lewis_2010_F16_lin_ex
%       (Single-design group) Implementation of the example in Sec. 5.2
%
%   vamvoudakis_lewis_2010_nonlin_ex
%       (Single-design group) Implementation of the example in Sec. 5.2
%
%   jiang_jiang_2014_engine_ex
%       (Single-design group) Implementation of the example in Sec. V. A.
%       *** NO LONGER SUPPORTED. This example is an unmatched uncertainty
%       example. Only the matched uncertainty algorithm alg_radp_matched.m
%       has been supported for this paper.
%
%   bian_jiang_2021_nonlin
%       (Single-design group) Implementation of the example in Sec. V. A.
%       *** NO LONGER SUPPORTED. This example uses a system for which the
%       control enters in the third power u^3, and a quartic control
%       penalty function. Not compatible with the code developed for the
%       paper.
%
% *************************************************************************

preset_group = 'CS_2ndorder_N1_3';
% preset_group = 'CS_2ndorder_N1_4';
% preset_group = 'CS_2ndorder_N1_7_no_x2p4';
% preset_group = 'CS_cip';

% preset_group = 'CS_2ndorder_N1_8';
% preset_group = 'CS_2ndorder_N1_7_no_x1p2';

% preset_group = 'vrabie_lewis_2009_hard_ex';
% preset_group = 'vamvoudakis_lewis_2010_F16_lin_ex';
% preset_group = 'vamvoudakis_lewis_2010_nonlin_ex';
% preset_group = 'jiang_jiang_2014_engine_ex';
% preset_group = 'bian_jiang_2021_nonlin';


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



