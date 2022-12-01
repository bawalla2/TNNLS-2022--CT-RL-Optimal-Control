function [figcount_out] = plotbode(sys_cell, wvec, plottypes, ttl_cell, lgd_text, axes_cell, figcount_in)

% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% plotbode.m
%
% Brent Wallace
%
% 1/23/20
%
% This function is meant to reduce the workload of plotting a large number
% of bode plots.
%
% CALL SYNTAX
%
% [figcount_out] = plotbode(sys_cell, wvec, plottypes, ttl, lgd_text, axes_cell, figcount_in)
%
% INPUTS
%
% sys_cell      -- A cell array containing the sytems to plot frequency
%               responses for. The systems will be plotted in the order of
%               their insertion in this cell array.
%               NOTE: If it is desired to plot multiple systems
%               corresponding to a list of a swept design parameter (e.g.,
%               plotting S and T on the same plot, sweeping w_g) then enter
%               this cell array as a cell array of cell arrays. E.g., if
%               S_cell and T_cell are, say, N-dimensional cell arrays
%               consisting of sensitivities and complementary sensitivities
%               (respectively), then to plot S and T on the same plot,
%               declare sys_cell = {S_cell , T_cell}.
% wvec          -- Can be either a vector or a cell array of vectors. If a 
%               vector is passed, all frequency responses will be evaluated 
%               at the frequencies in the vector. If a cell array is
%               passed, each frequency response will be evaluated the
%               frequencies of the vector in the corresponding cell array
%               entry
% plottypes     -- A 2-vector with entries [plotmag, plotph] where plotmag =
%               1 will plot magnitude responses, plotmag = 0 will not plot
%               magnitude responses. Likewise for plotph and phase
%               response.
% ttl_cell      -- A cell array containing titles for both plots. If the
%               empty cell is passed, then no title will be inserted. If a
%               1-element cell is passed, then the corresponding title will
%               be used for all plots.
% lgd_text      -- A cell array containing legend entries. If the empty
%               array is passed, then no legend will be inserted
% axes_cell     -- A cell array containing two 4-vectors containing axis window parameters. The
%               first vector corresponds to the axes parameters for the magnitude plot, the
%               second vector for the phase plot. If the
%               empty vector is passed in either of the two cells, then default axis windows will be
%               used for the respective plot. If a single 4-vector is
%               passed, then the same axes will be used for both plots.
% figcount_in   -- The current figure count at time of function call.
%
% OUTPUTS
%
% figount_out   -- The appropriate figure counter after function execution
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% *************************************************************************
%
% SETUP
%


% ***********************
%
% DETERMINE IF MULTIPLE SYSTEMS ARE TO BE PLOTTED FOR A SWEPT PARAMETER
%
% E.g., plotting S and T on the same plot, sweeping w_g. In this case, we
% say that we are "stacking" two systems into one.
%
if iscell(sys_cell{1})
   
    % If true, then muliple systems are to be plotted for the same sweep
    % parameter. The cell array already has the desired structure, so no
    % action is required.
    
    % Number of systems "stacked."
    stacknum = length(sys_cell);
    
else
    
    % Else, we are simply plotting for a list of systems. Embed this cell
    % array one layer deeper, so that the cell array has the required
    % structure.
    sys_cell = {sys_cell};
    
    % Number of systems "stacked" is 1.
    stacknum = 1;
    
end

% ***********************
%
% DETERMINE NUMBER OF SYSTEMS, NUMBER OF FREQUENCY POINTS
%

numsys = length(sys_cell{1});
numwpts = length(wvec);


% ***********************
%
% DETERMINE WHICH PLOTS ARE TO BE CREATED
%

plotmag = plottypes(1);     % 1 = plot magnitude response, 0 = don't plot
plotph = plottypes(2);      % 1 = plot phase response, 0 = don't plot


% ***********************
%
% DETERMINE IF MANUAL AXIS LIMITS ARE TO BE USED
%
if isempty(axes_cell)
    
    usemanaxes_mag = 0;        % If empty, use default axis windows
    usemanaxes_ph = 0;         % If empty, use default axis windows
    
elseif length(axes_cell) == 1   % Use same axes for both plots
    
    usemanaxes_mag = 1;
    usemanaxes_ph = 1;
    
    axes_mag = axes_cell{1};
    axes_ph = axes_cell{1};
    
else                            % Use axes specified in both cell elements
    
    usemanaxes_mag = 1;
    usemanaxes_ph = 1;
    
    axes_mag = axes_cell{1};
    axes_ph = axes_cell{2};
    
end

% ***********************
%
% DETERMINE TITLE OF PLOTS
%

if isempty(ttl_cell)            % Empty argument. No titles used.
    
    magttl = '';
    phttl = '';
    
elseif length(ttl_cell) == 1    % Use same title for both plots
        
    magttl = ttl_cell{1};
    phttl = ttl_cell{1};
    
else                            % Use both titles passed to function
    
    magttl = ttl_cell{1};
    phttl = ttl_cell{2};
    
end


% ***********************
%
% DETERMINE IF LEGEND IS TO BE USED
%

uselgd = ~isempty(lgd_text);        % If empty, do not use a legend

% ***********************
%
% DETERMINE FREQUENCIES OF EVALUATION
%

wvec_cell = cell(numsys,1);     % Cell array holding frequency vector of respective system

if ~iscell(wvec)    % Evaluate all systems at the same frequencies
    
    % Make sure the frequency vector is a column vector.
    wvec = wvec(:);
    
    for i = 1:numsys
       
        wvec_cell{i} = wvec;
        
    end
    
else                % Else, wvec is a cell array already
    
    % Make sure each vector in the cell array is a column vector.
    for i = 1:numsys
       
        tmp = wvec{i};
        tmp = tmp(:);       % Ensures wvec is a column vector.
        wvec_cell{i} = tmp;
        
    end
    
end

% ***********************
%
% INITIALIZE FIGURE COUNTER
%

figcount = figcount_in;

% *************************************************************************
%
% CALCULATE MAGNITUDE AND PHASE RESPONSES
%
% How this loop works:
%
% In the case that we are just plotting a list of systems, this loop simply
% runs through the systems and gets the frequency response data for each.
%
% In the case that we are plotting multiple "stacked" systems over N values
% of a sweep parameter (e.g., "stacking" S and T on the same plot, for N
% values of w_g) this loop calculates the frequency responses of each
% system in the stack first for the fixed parameter (e.g., it calculates
% frequency response data for S(1), then T(1)), and then runs the exterior
% loop down the direction of the sweep parameter (so, e.g., the order would
% be S(1), T(1), then S(2), T(2), ....). A single stack is treated as one
% legend entry item (so, e.g., the curves of S(1), T(1) will be treated as
% one curve and thus have one legend entry).
%

mag_cell = cell(numsys,1);        % Stores magnitude responses
ph_cell = cell(numsys,1);         % Stores phase responses

% Used for plotting only. The "effective" frequency vector that the
% frequency response data needs for plotting when multiple systems are
% stacked for a given parameter.
% E.g., if a list of systems is to be plotted, then wvec_eff_cell will be
% equal to wvec_cell. But if, e.g., S and T are to be put on the same plot,
% then the k-th entry of wvec_eff_cell will be 
%       [wvec_cell{k} ; nan ; wvec_cell{k}]
% This is still a vector, but it accommodates that the frequency response
% data has been stacked (with nan padding in between each system in the
% stack).
%
wvec_eff_cell = cell(numsys,1);


for syscount = 1:numsys

    for stackcount = 1:stacknum

        % Current cell array in the "stack."
        currcell = sys_cell{stackcount};

        % Current system in the current cell array.
        currsys = currcell{syscount};

        % Get frequency response data.
        [tmp_mag, tmp_ph] = bode(currsys,wvec_cell{syscount});

        % Vectorize the data.
        tmp_mag = squeeze(tmp_mag);
        tmp_ph = squeeze(tmp_ph);

        % Append this frequency response data to the current "stack" of
        % systems.
        if stackcount == 1

            mag = tmp_mag;
            ph = tmp_ph;
            wvec_eff = wvec_cell{syscount};

        else

            mag = [mag ; nan ; tmp_mag];
            ph = [ph ; nan ; tmp_ph];
            wvec_eff = [wvec_eff ; nan ; wvec_cell{syscount}];

        end

    end

    % Store frequency response data.
    mag_cell{syscount} = mag;
    ph_cell{syscount} = ph;
    wvec_eff_cell{syscount} = wvec_eff;

end


% *************************************************************************
%
% PLOT MAGNITUDE AND PHASE RESPONSES
%

% ***********************
%
% PLOT MAGNITUDE
%

if plotmag
   
    % ***********************
    %
    % PLOT
    %
    
    figure(figcount)

    for syscount = 1:numsys
        
        h = semilogx(wvec_eff_cell{syscount}, 20*log10(mag_cell{syscount})); 
        set(h, 'LineWidth', 2);
        
        if syscount == 1
            hold on
        end
        
    end
    
    
    % ***********************
    %
    % PLOT FORMATTING
    %
    
    % Set axes limits
    if usemanaxes_mag
        
        axis(axes_mag);
        
    end
    
    % Title
    title(magttl);
    
    % Axis labels
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude (dB)')
    
    % Grid
    grid on
    
    % Legend
    if uselgd
        
        lgd = legend(lgd_text);
        set(lgd, 'Location', 'Best');       % Put legend in empty spot
        
    end
    
    
    % Increment figure counter
    figcount = figcount + 1;
    
    
end



% ***********************
%
% PLOT PHASE
%

if plotph
   
    % ***********************
    %
    % PLOT
    %
    
    figure(figcount)
    
    for syscount = 1:numsys
        
        h = semilogx(wvec_eff_cell{syscount}, ph_cell{syscount});
        set(h, 'LineWidth', 2);
        
        if syscount == 1
            hold on
        end
        
    end
    
    
    % ***********************
    %
    % PLOT FORMATTING
    %
    
    % Set axes limits
    if usemanaxes_ph
        
        axis(axes_ph);
        
    end
    
    % Title
    title(phttl);
    
    % Axis labels
    xlabel('Frequency (rad/s)')
    ylabel('Phase (deg)')
    
    % Grid
    grid on
    
    % Legend
    if uselgd
        
        lgd = legend(lgd_text);
        set(lgd, 'Location', 'Best');       % Put legend in empty spot
        
    end
    
    
    % Increment figure counter
    figcount = figcount + 1;
    
    
end




% ***********************
%
% UPDATE FIGURE COUNTER OUTPUT
%

figcount_out = figcount;

end

