function figcount = plot_x_u(alg_settings_cell,...
                        out_data_cell, group_settings)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% PLOT STATE TRAJECTORY x(t) AND CONTROL SIGNAL u(t) 
%
% Brent Wallace  
%
% 2021-11-06
%
% This program, given a specified algorithm and data, plots state
% trajectory and control signal.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% plot_main(alg_settings_cell, out_data_cell, group_settings)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% alg_settings_cell     (Cell, each entry a Struct) Algorithm 
%                       settings/parameters for subsequent execution
%                       according to desired preset (fields vary by
%                       algorithm, see respective algorithm .m-file for
%                       details).
% out_data_cell         (Cell, each entry a Struct) Output data generated
%                       by the algorithm (fields vary by algorithm, see
%                       respective algorithm .m-file for details).
% group_settings        (Struct) contains preset group settings.
%                       Has the following fields:
%   savefigs            (Boolean) 1 = save figures to PDF. 0 = don't save.
%   relpath             (String) relative file path of folder to save plots
%                       to, if they are to be saved.
%   dolegend            (Boolean) 1 = include preset group legend on plots.
%                       0 = don't include legend.
%   lgnd                (Cell, each entry a string, optional) If dolegend =
%                       1, then this contains the legend to include in the
%                       plots.
%   sys_plot_settings   (Struct) contains system plot settings. See
%                       config_sys.m for details.
%   preset_group        (String) Tag of the current preset group being
%                       executed.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% 
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INITIALIZATION
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% Unpack plot settings
savefigs = group_settings.savefigs;
relpath = group_settings.relpath;
dolegend = group_settings.dolegend;
% do_individual_plots = group_settings.do_individual_plots;
sys_plot_settings = group_settings.sys_plot_settings;
figcount = group_settings.figcount;

% Number of designs to plot for
numpresets = size(alg_settings_cell,1);

% Extract system
sys = alg_settings_cell{1}.sys;

% Get total system order. This includes the n-state variables, as well as
% any unknown dynamic variables
sys_order = size(out_data_cell{1}.xmat, 2);

% State trajectory x(t) unit scaling. E.g., if x_1(t) is angular
% displacement but is desired in degrees, declare the field 'sclvec' and
% put sclvec(1) = 180/pi.
x_sclvec = sys_plot_settings.x_sclvec;

% Time axis label
tlabel = group_settings.tlabel;

% Legend for all states in one plot
x_t_state_lgd = sys_plot_settings.x_t_state_lgd;

% ***********************
%
% EXTRACT LEGEND IF USER SPECIFIED
%  

% Extract legend, if user specified
if dolegend
    lgnd = group_settings.lgnd;
end


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
% STATE TRAJECTORY PLOTS
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
% *************************************************************************
%
% INDIVIDUAL STATES
% 
% *************************************************************************
% *************************************************************************

% Keeps track of how to window the time axis
tmax = -1;

for i = 1:sys_order
   
    % Open new figure
    figure(figcount);
    
    % Preliminary settings
    hold on
    
    % Extract scaling for x_i
    x_scl = x_sclvec(i);
    
    % Plot state trajectory x_i(t) for each preset in the group
    for j = 1:numpresets
        
        % Extract data from current preset
        tvec = out_data_cell{j}.tvec;
        x_i_t = out_data_cell{j}.xmat(:,i);
        
        % Plot
        h_fig = plot(tvec, x_scl * x_i_t);
        
        % Keep track of how to window time axis
        tmax = max(tmax, tvec(end));
        
    end
    
    % Title, labels
    title(sys_plot_settings.x_t_title{i});
    xlabel(tlabel);
    xlim([0 tmax]);
    ylabel(sys_plot_settings.x_t_ylabel{i});
    
    % Legend
    if dolegend
        lgd = legend(lgnd);             % Create legend
    end
       
    % Format plot
    p_sett.figcount = figcount;
    plot_format(p_sett, group_settings); 
    
    % SAVE PLOT
    if savefigs
        filename = sys_plot_settings.x_t_filename{i};
        savepdf(figcount, relpath, filename); 
    end
    
    % Increment figure counter
    figcount = figcount + 1; 
    
end

% *************************************************************************
% *************************************************************************
%
% ALL STATES
% 
% *************************************************************************
% *************************************************************************


% Open new figure
figure(figcount);

% Preliminary settings
hold on

% Plot state trajectory x(t) for each preset in the group
for j = 1:numpresets

    % Extract data from current preset
    tvec = out_data_cell{j}.tvec;
    x_t = out_data_cell{j}.xmat;
    
    % Scale 
    for i = 1:sys.n
        x_t(:,i) = x_t(:,i) * x_sclvec(i);
    end

    % Plot
    h_fig = plot(tvec, x_t);

end

% Title, labels
title('State Trajectory $x(t)$');
xlabel(tlabel);
xlim([0 tmax]);
ylabel('$x(t)$');

% Legend
lgd = legend(x_t_state_lgd);

% Format plot
p_sett.figcount = figcount;
plot_format(p_sett, group_settings); 

% SAVE PLOT
if savefigs
    filename = 'x_t_all';
    savepdf(figcount, relpath, filename); 
end

% Increment figure counter
figcount = figcount + 1; 



%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CONTROL SIGNAL PLOTS
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************


for i = 1:sys.m
   
    % Open new figure
    figure(figcount);
    
    % Preliminary settings
    hold on
    
    % Plot control u(t) for each preset in the group
    for j = 1:numpresets
        
        % Extract data from current preset
        tvec = out_data_cell{j}.tvec;
        u_i_t = out_data_cell{j}.umat(:,i);
        
        % Plot
        h_fig = plot(tvec, u_i_t);

    end
    
    % Title, axes labels
    title(sys_plot_settings.u_t_title{i});
    xlabel(tlabel);
    xlim([0 tmax]);
    ylabel(sys_plot_settings.u_t_ylabel{i});
    
    % Legend
    if dolegend
        lgd = legend(lgnd);             % Create legend
    end
    
    % Format plot
    p_sett.figcount = figcount;
    plot_format(p_sett, group_settings); 
    
    % SAVE PLOT
    if savefigs
        filename = sys_plot_settings.u_t_filename{i};
        savepdf(figcount, relpath, filename); 
    end
    
    % Increment figure counter
    figcount = figcount + 1; 
    
end

