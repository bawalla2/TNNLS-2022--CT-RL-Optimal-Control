% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% MERGE PRESET GROUP SETTINGS/DATA FROM TWO PREVIOUSLY EXECUTED x_0 SWEEP
% PRESET GROUPS
%
% Brent Wallace  
%
% 2022-03-25
%
% This program, given relative paths to two previously generated x_0 sweep
% preset group 'alg_settings_cell', 'out_data_cell', and 'group_settings'
% objects, merges them into one new set of 'alg_settings_cell',
% 'out_data_cell', and 'group_settings' objects.
%
% NOTE: The two previously executed x_0 sweep preset groups must have been
% executed for the same x_0 sweep values. In addition, this code only
% supports at maximum one set of design parameters per algorithm. E.g., if
% the first preset group executed IRL, and the second group exectued IRL,
% the merged group could only contain the IRL presets from at most one of
% the parent groups.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% merge_x0_sweeps
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
% Settings are changed manually in this script. However, the
% following variables will need to be declared manually to run the script:
%
%
% relpath_out           (String) Relative path to write the merged
%                       'alg_settings_cell', 'out_data_cell', and
%                       'group_settings' structs to.
%
% inherit_1_1_2_0       (Bool) Inherit the default properties of the
%                       first 'group_settings' object (=1) (i.e., the one
%                       pointed to by 'relpath_1'), or from the second
%                       group_settings' object (=0).
%
% relpath_1, relpath_2    (Strings) Relative paths to the
%                       'alg_settings_cell', 'out_data_cell', and
%                       'group_settings' structs corresponding to the first
%                       and second x_0 sweeps, respectively.
%
% alg_list_1, alg_list_2  (Cell arrays, each entry a string) An array
%                       consisting of the algorithms from each x_0 sweep
%                       for which it is desired to merge into the output
%                       preset data. See above note; these lists must be
%                       disjoint. The algorithm strings can take on the
%                       following values:
%   'irl'               IRL
%   'spi'               SPI
%   'radp_matched'      RADP (matched uncertainty)
%   'vi'                VI
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% This function outputs the following, which are written to the relative
% path determined by 'relpath_out':
%
% alg_settings_cell     (Cell array, each entry a struct) The merged
%                       'alg_settings_cell' object from the two input
%                       'alg_settings_cell' objects.
% out_data_cell         (Cell array, each entry a struct) The merged
%                       'out_data_cell' object from the two input
%                       'out_data_cell' objects.
% group_settings        (Struct) The merged 'group_settings' object from
%                       the two input 'group_settings' objects. 
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
% RELATIVE PATHS TO WRITE MERGED DATA TO 
% 
% *************************************************************************

% Write to the figures output folder
relpath_out = '00 figures/';

% Add a timestamp
time = fix(clock);
timestamp = '';
for i = 1:length(time)
   timestamp = [timestamp , num2str(time(i))];
   if i < length(time)
       timestamp = [timestamp , '_'];
   end
end
timestamp = [timestamp, '\'];
relpath_out = [relpath_out, timestamp];   % Update relative path

% *************************************************************************
%
% WHICH PRESET GROUP TO INHERIT DEFAULT PROPERTIES FROM
%
% First preset group:       = 1
% Second preset group:      = 0
% 
% *************************************************************************

inherit_1_1_2_0 = 1;    % Inherit from first preset group
% inherit_1_1_2_0 = 0;    % Inherit from second preset group

% *************************************************************************
%
% RELATIVE PATHS TO PRE-EXISTING DATA
% 
% *************************************************************************

% ***********************
%       
% RELATIVE PATH TO DATASET 1
%    

% % CS I -- V1 -- e_1(t)
% relpath_1 = '01 data/CS_2ndorder_N1_3/V1/e_1_t/';
% % CS I -- V2 -- e_1(t)
% relpath_1 = '01 data/CS_2ndorder_N1_3/V2/e_1_t/';
% % CS I -- V1 -- e_2(t)
% relpath_1 = '01 data/CS_2ndorder_N1_3/V1/e_2_t/';
% % CS I -- V1 -- e_3(t)
% relpath_1 = '01 data/CS_2ndorder_N1_3/V1/e_3_t/';

% % CS II -- V1 -- e_1(t)
% relpath_1 = '01 data/CS_2ndorder_N1_4/V1/e_1_t/';
% % CS II -- V1 -- e_2(t)
% relpath_1 = '01 data/CS_2ndorder_N1_4/V1/e_2_t/';
% % CS II -- V1 -- e_3(t)
% relpath_1 = '01 data/CS_2ndorder_N1_4/V1/e_3_t/';

% CS III -- V2 -- IRL, RADP
relpath_1 = '01 data/CS_2ndorder_N1_7_no_x2p4/V2/IRL_RADP/';


% ***********************
%       
% RELATIVE PATH TO DATASET 2
%    

% % CS I -- V1 -- e_1(t)
% relpath_2 = '01 data/CS_2ndorder_N1_3/V1/e_1_t/';
% % CS I -- V2 -- e_1(t)
% relpath_2 = '01 data/CS_2ndorder_N1_3/V2/e_1_t/';
% % CS I -- V2 -- e_2(t)
% relpath_2 = '01 data/CS_2ndorder_N1_3/V2/e_2_t/';
% % CS I -- V2 -- e_3(t)
% relpath_2 = '01 data/CS_2ndorder_N1_3/V2/e_3_t/';


% % CS II -- V2 -- e_1(t)
% relpath_2 = '01 data/CS_2ndorder_N1_4/V2/e_1_t/';
% % CS II -- V2 -- e_2(t)
% relpath_2 = '01 data/CS_2ndorder_N1_4/V2/e_2_t/';
% % CS II -- V2 -- e_3(t)
% relpath_2 = '01 data/CS_2ndorder_N1_4/V2/e_3_t/';

% CS III -- V2 -- IRL, RADP, VI
relpath_2 = '01 data/CS_2ndorder_N1_7_no_x2p4/V2/IRL_RADP_VI/';


% *************************************************************************
%
% ALGORITHM LISTS OF EACH DATASET TO MERGE
% 
% *************************************************************************


% ***********************
%       
% ALGORITHM LIST OF DATASET 1
%  

alg_list_1 =  {
                'irl'
%                 'spi'
                'radp_matched'
%                 'vi'
                               };
                           
                           
% ***********************
%       
% ALGORITHM LIST OF DATASET 2
%  

alg_list_2 =  {
%                 'irl'
%                 'spi'
%                 'radp_matched'
                'vi'
                               };                           

                           
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
% EXTRACT DATA, PREPARE FOR MERGING
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% *************************************************************************
%
% EXTRACT PRE-EXISTING DATA
% 
% *************************************************************************

% ***********************
%       
% LOAD DATASET 1
% 

% % Load data -- alg_settings_cell struct 1
% data = load([relpath_1 'alg_settings_cell.mat']);
% alg_settings_cell_1 = data.alg_settings_cell;
% 
% % Load data -- out_data_cell struct 1
% data = load([relpath_1 'out_data_cell.mat']);
% out_data_cell_1 = data.out_data_cell;

% Load data -- group_settings struct 1
data = load([relpath_1 'group_settings.mat']);
group_settings_1 = data.group_settings;

% ***********************
%       
% LOAD DATASET 2
% 

% % Load data -- alg_settings_cell struct 2
% data = load([relpath_2 'alg_settings_cell.mat']);
% alg_settings_cell_2 = data.alg_settings_cell;
% 
% % Load data -- out_data_cell struct 2
% data = load([relpath_2 'out_data_cell.mat']);
% out_data_cell_2 = data.out_data_cell;

% Load data -- group_settings struct 1
data = load([relpath_2 'group_settings.mat']);
group_settings_2 = data.group_settings;

% *************************************************************************
%
% GET PRESET GROUP SETTINGS
% 
% *************************************************************************

% Number of ICs in the sweep
numICs = group_settings_1.numICs;

% *************************************************************************
%
% GET INDICES CORRESPONDING TO THE SPECIFIED ALGORITHMS FROM EACH PRESET
% 
% *************************************************************************

for i = 1:2
   
    switch i
        % Find indices of preset group 1
        case 1
            group_settings_i = group_settings_1;
            alg_list_i = alg_list_1;
        % Find indices of preset group 2    
        case 2
            group_settings_i = group_settings_2;
            alg_list_i = alg_list_2;
    end
    
    % Indices of presets corresponding to specified algorithms.
    inds_i = [];
    
    for j = 1:size(alg_list_i, 1)
        
        % Get current algorithm
        alg = alg_list_i{j};
        
        % Find this algorithm's index in the preset group
        ind_alg = -1;        
        for k = 1:group_settings_i.numalgs
            if strcmp(alg, group_settings_i.alg_list{k})
                ind_alg = k;
            end
        end
        
        % If algorithm not found, throw an error
        if ind_alg == -1
           error([' *** ERROR: ALGORITHM    ' alg '     NOT FOUND '...
               'IN PRESET GROUP   ' num2str(i) ])
        end
        
        % Add the indices in the preset group corresponding to this
        % algorithm
        inds_i =    [   inds_i
                        (ind_alg-1) * numICs + (1:numICs)'  ];
        
    end
    
    % Store the indices
    switch i
        case 1
            inds_1 = inds_i; 
        case 2
            inds_2 = inds_i;
    end    
    
end


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% MERGE DATA
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% *************************************************************************
%
% MERGE 'alg_settings_cell' and 'out_data_cell' OBJECTS
% 
% *************************************************************************

[alg_settings_cell, out_data_cell] = merge_settings_data(relpath_1,...
                                        inds_1, relpath_2, inds_2);
                                    
                                    
% *************************************************************************
%
% MERGE 'group_settings' OBJECTS
% 
% *************************************************************************     

if inherit_1_1_2_0

    % Inherit from the first 'group_settings' object
    group_settings = group_settings_1;

else
    
    % Inherit from the second 'group_settings' object
    group_settings = group_settings_2;

end
                                    
% ***********************
%       
% FIELDS TO MODIFY IN THE MERGED 'group_settings' OBJECT TO REFLECT THE
% MERGE
% 

% Algorithm list
alg_list_merged = [     alg_list_1
                        alg_list_2  ];

% Number of algorithms in the sweep
numalgs_merged = size(alg_list_merged, 1);
                    
% Number of presets in the preset group    
numpresets_merged = numICs * numalgs_merged;

% ***********************
%       
% MODIFY MERGED 'group_settings' OBJECT TO REFLECT THE MERGE
% 

% Algorithm list
group_settings.alg_list = alg_list_merged;

% Number of algorithms in the sweep
group_settings.numalgs = numalgs_merged;
                    
% % Number of presets in the preset group    
% group_settings.numpresets = numpresets_merged;


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% WRITE MERGED DATA
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% Make directory to save data to
mkdir(relpath_out);              

% Save data -- alg_settings_cell struct
varname = 'alg_settings_cell';
save([relpath_out varname], varname);

% Save data -- out_data_cell struct
varname = 'out_data_cell';
save([relpath_out varname], varname);

% Save data -- group_settings struct
varname = 'group_settings';
save([relpath_out varname], varname);
