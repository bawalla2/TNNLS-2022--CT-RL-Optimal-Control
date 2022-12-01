function plot_format(plot_settings, group_settings)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% APPLY BASIC FORMATTING TO PLOT
%
% Brent Wallace  
%
% 2022-01-13
%
% This program applies basic formatting to a given plot.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% plot_format(figcount, group_settings)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% plot_settings     (Struct) Structure containing various plot controls.
%                   Has the following fields:
%   figcount        (Integer) Figure number of the figure to edit
%                   formatting of.
%
% group_settings    (Struct) Structure containing the overall preset group
%                   settings. See main.m for description of fields.
%
% ***** TO OVERRIDE DEFAULT SETTINGS
%
% Each of the formatting options available in this program can be found
% under the 'settings_cell' variable (declared below). If it is desired to
% override any of these formatting values, then in the 'plot_settings'
% input struct, declare a field named '.custom_sett' under which to
% override them. To illustrate, here is an example:
%
%   If it is desired to change the default title font size, the
%   corresponding formatting variable is 'ttl_fontsize'. If one wants to
%   change this to, say, 10, they would in plot_settings.custom_sett
%   declare:
%               plot_settings.custom_sett.ttl_fontsize = 10;
%
% NOTE: Do not declare '.custom_sett' if no default settings are to be
% overriden.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% NONE (Re-formatted plots)
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INITIALIZATION
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
%
% UNPACK USER SETTINGS
% 
% *************************************************************************

% ***********************
%       
% GET FIGURE SPECIFIED BY USER
%  

% Figure number of figure to be edited
figcount = plot_settings.figcount;

% Get figure
fig = figure(figcount);

% Get CurrentAxes object of figure
ax = fig.CurrentAxes;

% Get x, y, z axis objects
ax_x = ax.XAxis;
ax_y = ax.YAxis;
ax_z = ax.ZAxis;

% Get any line object present in the figure
lines = findobj(fig, 'Type', 'Line');
haslines = ~isempty(lines);

% Get Legend object of figure (is inside the CurrentAxes object)
lgd = ax.Legend;
haslgd = ~isempty(lgd);

% Get if the plot is a 2D plot or 3D plot
twoD1_threeD0 = isequal(ax.View, [0 90]);

% Get x, y, z window limits
ax_x_lim = ax_x.Limits;
ax_y_lim = ax_y.Limits;
ax_z_lim = ax_z.Limits;

% Get x, y, z axis ticks
ax_x_ticks = ax_x.TickValues;
ax_y_ticks = ax_y.TickValues;
ax_z_ticks = ax_z.TickValues;

% ***********************
%       
% CHECK FOR CUSTOM SETTINGS SPECIFIED BY USER
%  

% Any custom settings that have been specified by the user to override the
% defaults
override_sett = isfield(plot_settings, 'custom_sett');
if override_sett   
    custom_sett = plot_settings.custom_sett;
end

% *************************************************************************
%
% CELL OF ALL FORMATTING SETTING VARIABLE NAMES
% 
% *************************************************************************

settings_cell = {
                    % Title formatting
                    'ttl_font'
                    'ttl_fontsize'
                    'ttl_fontweight'
                    'ttl_interp'
                    % Axis tick formatting
                    'tick_font'
                    'tick_fontsize'
                    'tick_fontweight'
                    'tick_interp'  
                    'tick_rot'
                    'keep_ticks'
                    % Axis label formatting
                    'lbl_font'
                    'lbl_fontsize'
                    'lbl_fontweight'
                    'lbl_interp' 
                    % Legend formatting
                    'lgd_font'
                    'lgd_fontsize'
                    'lgd_fontweight'
                    'lgd_interp'   
                    'lgd_loc'
                    % Line formatting
                    'line_width'
                    % Grid formatting
                    'grid_on'
                                    };
                                
len_settings_cell = size(settings_cell, 1);                                

% *************************************************************************
%
% DEFAULT PLOT FORMATTING SETTINGS
% 
% *************************************************************************

% ***********************
%       
% TITLE FORMATTING
%  

% Title font -- MATLAB Default: 'Helvetica'
ttl_font = 'Helvetica';

% Title font size -- MATLAB Default: 11
ttl_fontsize = 15;

% Title font weight -- MATLAB Default: 'bold'
ttl_fontweight = 'bold';

% Title interpreter -- MATLAB Default: 'tex'
ttl_interp = 'latex';


% ***********************
%       
% AXIS TICK FORMATTING
%  

% Tick font -- MATLAB Default: 'Helvetica'
tick_font = 'Helvetica';

% Tick font size -- MATLAB Default: 10
tick_fontsize = 12;

% Tick font weight -- MATLAB Default: 'normal'
tick_fontweight = 'bold';

% Tick interpreter -- MATLAB Default: 'tex'
tick_interp = 'tex';

% Tick label rotation -- MATLAB Default: 0 (vertical)
tick_rot = 0;

% Keep current tick values (=1), or let formatting changes dictate new tick
% values (=0)
keep_ticks = 1;


% ***********************
%       
% AXIS LABEL FORMATTING
%  

% Label font -- MATLAB Default: 'Helvetica'
lbl_font = 'Helvetica';

% Label font size -- MATLAB Default: 11
lbl_fontsize = 20;

% Label font weight -- MATLAB Default: 'normal'
lbl_fontweight = 'bold';

% Label interpreter -- MATLAB Default: 'tex'
lbl_interp = 'latex';


% ***********************
%       
% LEGEND FORMATTING
%  

% Legend font -- MATLAB Default: 'Helvetica'
lgd_font = 'Helvetica';

% Legend font size -- MATLAB Default: 9
lgd_fontsize = 12;

% Legend font weight -- MATLAB Default: 'normal'
lgd_fontweight = 'normal';
% lgd_fontweight = 'bold';

% Legend interpreter -- MATLAB Default: 'tex'
lgd_interp = 'latex';

% Legend location -- MATLAB Default: 'northeast'
if twoD1_threeD0
    lgd_loc = 'northeast';
else
    lgd_loc = 'northeast';
end


% ***********************
%       
% LINE FORMATTING
%

% Line width -- MATLAB Default: 0.5
line_width = 2;

% ***********************
%       
% GRID FORMATTING
% 

% Turn grid on or not -- MATLAB Default: 0
grid_on = 1;


% *************************************************************************
%
% OVERRIDE DEFAULT SETTINGS IF USER SPECIFIES
% 
% *************************************************************************

if override_sett

    for i = 1:len_settings_cell
        
        currsett = settings_cell{i};
        
        if isfield(plot_settings.custom_sett, currsett)
            
            switch currsett

                
                % ***********************
                %       
                % TITLE FORMATTING
                % 
                
                case 'ttl_font'
                    
                    ttl_font = custom_sett.ttl_font;                
                
                case 'ttl_fontsize'
                    
                    ttl_fontsize = custom_sett.ttl_fontsize;
                    
                case 'ttl_fontweight'
                    
                    ttl_fontsize = custom_sett.ttl_fontweight;
                        
                case 'ttl_interp'
                    
                    ttl_interp = custom_sett.ttl_interp;   

                % ***********************
                %       
                % AXIS TICK FORMATTING
                %                     
                    
                case 'tick_font'
                    
                    tick_font = custom_sett.tick_font;  
                
                case 'tick_fontsize'
                    
                    tick_fontsize = custom_sett.tick_fontsize;
                    
                case 'tick_fontweight'
                    
                    tick_fontsize = custom_sett.tick_fontweight;
                        
                case 'tick_interp'
                    
                    tick_interp = custom_sett.tick_interp;        
                    
                case 'tick_rot'
                    
                    tick_rot = custom_sett.tick_rot;                     

                case 'keep_ticks'
                    
                    keep_ticks = custom_sett.keep_ticks;                      
                    
                % ***********************
                %       
                % AXIS LABEL FORMATTING
                % 
                
                case 'lbl_font'
                    
                    lbl_font = custom_sett.lbl_font;  
                
                case 'lbl_fontsize'
                    
                    lbl_fontsize = custom_sett.lbl_fontsize;
                    
                case 'lbl_fontweight'
                    
                    lbl_fontsize = custom_sett.lbl_fontweight;
                        
                case 'lbl_interp'
                    
                    lbl_interp = custom_sett.lbl_interp;                       
                    
                % ***********************
                %       
                % LEGEND FORMATTING
                %                     
                    
                case 'lgd_font'
                    
                    lgd_font = custom_sett.lgd_font;                
                
                case 'lgd_fontsize'
                    
                    lgd_fontsize = custom_sett.lgd_fontsize;
                    
                case 'lgd_fontweight'
                    
                    lgd_fontsize = custom_sett.lgd_fontweight;
                        
                case 'lgd_interp'
                    
                    lgd_interp = custom_sett.lgd_interp;
                    
                case 'lgd_loc'
                    
                    lgd_loc = custom_sett.lgd_loc;    
                    
                % ***********************
                %       
                % LINE FORMATTING
                %                    

                case 'line_width'
                    
                    line_width = custom_sett.line_width;                 
                
                % ***********************
                %       
                % GRID FORMATTING
                %                     
                    
                case 'grid_on'
                    
                    grid_on = custom_sett.grid_on;                       
                    
                % ***********************
                %
                % THROW ERROR IF TAG DOES NOT COME UP A MATCH
                %   

                otherwise

                    error(['*** ERROR: PLOT SETTING'...
                        ' VARIABLE NAME NOT RECOGNIZED ***']);  
                    
                    
            end
            
        end
        
    end
    
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



% ***********************
%       
% TITLE FORMATTING
%  

% Title font
set(ax.Title, 'FontName', ttl_font);

% Title font size
set(ax.Title, 'FontSize', ttl_fontsize);

% Title font weight
set(ax.Title, 'FontWeight', ttl_fontweight);

% Title interpreter
set(ax.Title, 'interpreter', ttl_interp);


% ***********************
%       
% AXIS TICK FORMATTING
%  

% Tick font
set(ax_x, 'FontName', tick_font);
set(ax_y, 'FontName', tick_font);
set(ax_z, 'FontName', tick_font);

% Tick font size
set(ax_x, 'FontSize', tick_fontsize);
set(ax_y, 'FontSize', tick_fontsize);
set(ax_z, 'FontSize', tick_fontsize);

% Tick font weight
set(ax_x, 'FontWeight', tick_fontweight);
set(ax_y, 'FontWeight', tick_fontweight);
set(ax_z, 'FontWeight', tick_fontweight);

% Tick interpreter
set(ax_x, 'TickLabelInterpreter', tick_interp);
set(ax_y, 'TickLabelInterpreter', tick_interp);
set(ax_z, 'TickLabelInterpreter', tick_interp);

% Tick label rotation
set(ax_x, 'TickLabelRotation', tick_rot);
set(ax_y, 'TickLabelRotation', tick_rot);
set(ax_z, 'TickLabelRotation', tick_rot);

% Return to original axes limits
set(ax_x, 'Limits', ax_x_lim);
set(ax_y, 'Limits', ax_y_lim);
set(ax_z, 'Limits', ax_z_lim);

% Keep current tick values (=1), or let formatting changes dictate new tick
% values (=0)
if keep_ticks
    set(ax_x, 'TickValues', ax_x_ticks);
    set(ax_y, 'TickValues', ax_y_ticks);
    set(ax_z, 'TickValues', ax_z_ticks);
end


% ***********************
%       
% AXIS LABEL FORMATTING
%  

% Label font
set(ax.XLabel, 'FontName', lbl_font);
set(ax.YLabel, 'FontName', lbl_font);
set(ax.ZLabel, 'FontName', lbl_font);

% Label font size
set(ax.XLabel, 'FontSize', lbl_fontsize);
set(ax.YLabel, 'FontSize', lbl_fontsize);
set(ax.ZLabel, 'FontSize', lbl_fontsize);

% Label font weight
set(ax.XLabel, 'FontWeight', lbl_fontweight);
set(ax.YLabel, 'FontWeight', lbl_fontweight);
set(ax.ZLabel, 'FontWeight', lbl_fontweight);

% Label interpreter
set(ax.XLabel, 'interpreter', lbl_interp);
set(ax.YLabel, 'interpreter', lbl_interp);
set(ax.ZLabel, 'interpreter', lbl_interp);


% ***********************
%       
% LEGEND FORMATTING
%  

if haslgd

    % Legend font
    set(lgd, 'FontName', lgd_font);

    % Legend font size
    set(lgd, 'FontSize', lgd_fontsize);

    % Legend font weight
    set(lgd, 'FontWeight', lgd_fontweight);

    % Legend interpreter
    set(lgd, 'interpreter', lgd_interp);

    % Legend location
    set(lgd, 'location', lgd_loc);
    

end



% ***********************
%       
% LINE FORMATTING
%  

if haslines
   
    % Line width
    set(lines(:), 'LineWidth', line_width);
    
end

% ***********************
%       
% GRID FORMATTING
%  

% Turn grid on or not
if grid_on
    grid on;
end


