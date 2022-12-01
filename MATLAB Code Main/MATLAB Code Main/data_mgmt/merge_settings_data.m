function [alg_settings_cell, out_data_cell] = merge_settings_data(...
                relpath_1, inds_1, relpath_2, inds_2)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% MERGE PRESET GROUP SETTINGS/DATA FROM TWO PREVIOUSLY EXECUTED PRESET
% GROUPS
%
% Brent Wallace  
%
% 2022-03-25
%
% This program, given relative paths to two previously generated preset
% group settings and output data cells, merges them to create a new
% 'alg_settings_cell' and 'out_data_cell' with entries from the indices
% specified in the index vectors 'inds_1', 'inds_2'.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% [alg_settings_cell, out_data_cell] = merge_settings_data(relpath_1,...
%                                         inds_1, relpath_2, inds_2);
%
% NOTE: Before running this program, make sure the MATLAB directory is at
% the folder housing 'main.m'
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% relpath_1              (String) Relative path to the first
%                       'alg_settings_cell' and 'out_data_cell' structs to
%                       merge.
% inds_1                 (Vector of integers) The indices of the first
%                       preset group to collect for the merge.
% relpath_2              (String) Relative path to the second
%                       'alg_settings_cell' and 'out_data_cell' structs to
%                       merge.
% inds_2                 (Vector of integers) The indices of the second
%                       preset group to collect for the merge.
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
% EXTRACT PRE-EXISTING DATA
% 
% *************************************************************************

% Load data -- alg_settings_cell struct 1
data = load([relpath_1 'alg_settings_cell.mat']);
alg_settings_cell_1 = data.alg_settings_cell;

% Load data -- out_data_cell struct 1
data = load([relpath_1 'out_data_cell.mat']);
out_data_cell_1 = data.out_data_cell;

% Load data -- alg_settings_cell struct 2
data = load([relpath_2 'alg_settings_cell.mat']);
alg_settings_cell_2 = data.alg_settings_cell;

% Load data -- out_data_cell struct 2
data = load([relpath_2 'out_data_cell.mat']);
out_data_cell_2 = data.out_data_cell;


% *************************************************************************
%
% FILTER PRE-EXISTING DATA BY SPECIFIED INDICES
% 
% *************************************************************************

% Filter data -- alg_settings_cell struct 1
alg_settings_cell_1 = alg_settings_cell_1(inds_1);

% Filter data -- out_data_cell struct 1
out_data_cell_1 = out_data_cell_1(inds_1);

% Filter data -- alg_settings_cell struct 2
alg_settings_cell_2 = alg_settings_cell_2(inds_2);

% Filter data -- out_data_cell struct 1
out_data_cell_2 = out_data_cell_2(inds_2);


% *************************************************************************
%
% MERGE DATA INTO NEW CELL ARRAYS
% 
% *************************************************************************

% Merge data -- alg_settings_cell
alg_settings_cell = [   alg_settings_cell_1
                        alg_settings_cell_2     ];
                    
% Merge data -- out_data_cell
out_data_cell = [   out_data_cell_1
                    out_data_cell_2     ];
