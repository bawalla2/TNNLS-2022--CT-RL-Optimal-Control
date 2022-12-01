function nt = eval_noise(t, m, noise)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% EVALUATE PROBING NOISE n(t)
%
% Brent Wallace  
%
% 2021-11-06
%
% This program, given a tag of a desired noise signal, input dimension m,
% and time t evaluates the probing noise n(t).
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% nt = eval_noise(t, m, noise)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
%
% t         (Double) Time (sec)
% m         (Integer) System input dimension
% noise     (Struct) Contains properties of the pecific noise signal to be
%           evaluated. Has the following fields:
%   tag     (String) Tag of the noise signal to be evaluated
%
% NOTE: If the noise is of tag 'sum_sinusoids', additional fields are
% required (see description above the respective tag).
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% nt        (m-dim. vector) Evaluation of probing noise n(t)
%
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% BEGIN MAIN
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

switch noise.tag

    % ***********************
    %
    % ZERO NOISE SIGNAL
    %
    
    case '0'
        
        nt = zeros(m,1);
        
    % ***********************
    %
    % SUM OF SINES AND COSINES
    %
    % Requires the fields in the noise struct:
    %
    %   cos1_sin0       M-dim. binary vector, where M denotes the total
    %                   number of sinusoids in the noise signal. For the
    %                   ith sinusoid (1 <= i <= M) 1 = make cos(), 0 = make
    %                   sin().
    %   scalevec        M-dim. vector. For ith sinusoid, is the scalar
    %                   multiple to put in front of the ith sinusoid
    %   freqvec         M-dim. vector. For ith sinusoid, is the frequency
    %                   of the ith sinusoid
    %   
    % E.g., to make the noise n(t) = cos(t) - 2*sin(5t), we would have
    % cos1_sin0 = [1;0], scalevec = [1;-2], freqvec = [1;5].
    %   
    %   
    
    case 'sum_sinusoids'
        
        M = size(noise.cos1_sin0, 1);   % Number of sinusoids in the signal
        nt = zeros(m,1);
        
        for i = 1:M
            a = noise.scalevec(i);
            w = noise.freqvec(i);
            if noise.cos1_sin0(i)
                nt = nt + a * ones(m,1) * cos(w * t);
            else    
                nt = nt + a * ones(m,1) * sin(w * t);
            end 
        end

    % ***********************
    %
    % DECAYING SINUSOIDAL SIGNAL
    %
    
    case 'dacaying_sinusoid'   
        
        a = 5;              % Amplitude
        w = 1;              % Frequency
        alpha = 0.01;       % Decay rate
        
        nt = ones(m, 1) * a * cos(w * t) * exp(-alpha * t);
        
    % ***********************
    %
    % K.G. VAMVOUDAKIS, F.L. LEWIS (2010) -- F16 LINEAR EXAMPLE
    %
    
    case 'vamvoudakis_lewis_2010_ex'
        
        nt = ones(m,1) * exp(-0.001*t)*37*(sin(t)^2*cos(t) +...
            sin(2*t)^2*cos(0.1*t) +...
            sin(-1.2*t)^2*cos(0.5*t) + sin(t)^5 + sin(1.12*t)^2 +...
            cos(2.4*t)*sin(2.4*t)^3);
 
    % ***********************
    %
    % K.G. VAMVOUDAKIS, F.L. LEWIS (2010) -- 2ND ORDER NONLINEAR
    %
    % NOTE: No information is available on the actual probing noise used.
    % This noise is derived from the F16 noise provided and has been
    % modified to recreate the results of the paper.
    %
    
    case 'vamvoudakis_lewis_2010_ex_nonlin'
        
        nt = ones(m,1) * exp(-0.001*t)*1*(sin(t)^2*cos(t) +...
            sin(2*t)^2*cos(0.1*t) +...
            sin(-1.2*t)^2*cos(0.5*t) + sin(t)^5 + sin(1.12*t)^2 +...
            cos(2.4*t)*sin(2.4*t)^3);        
 
        
    % ***********************
    %
    % Y. JIANG, Z.-P. JIANG (2014) -- JET ENGINE EXAMPLE
    %
    
    case 'jiang_jiang_2014_engine_ex'
        
%         % As shown in paper
%         nt = ones(m, 1) * 10 * cos(0.1 * t);   
        
        % As implemented in code
        nt = ones(m, 1) * (-10 * cos(0.1 * t) + 1 * cos(t));    
        
%         nt = ones(m, 1) * (-10 * cos(0.1 * t) + 1 * cos(t) + 3 * cos(2*t)); 
%         nt = ones(m, 1) * 1 * cos(0.1 * t); 
%         nt = ones(m, 1) * 10 * cos(0.1 * t) * exp(-0.1 * t); 
%         nt = 0;
        
    % ***********************
    %
    % T. BIAN, Z.-P. JIANG (2021) -- NONLINEAR EXAMPLE
    %
    
    case 'bian_jiang_2021_nonlin'
        
        % As originally in code
        nt = ones(m, 1) * 1 * sin(5 * t);
        
%         % Modified
% %         nt = ones(m, 1) * ((1)*sin(7*t) + (-0.5)*cos(10*t));
%         nt = ones(m, 1) * ((1)*sin(7*t) + (-1)*cos(10*t));

    % ***********************
    %
    % RADP ON 2ND ORDER SYSTEM -- SWEEP x_0
    %   
    
    case 'radp_2ndorder_sweep_x0'
        
        nt = ones(m, 1) * (20) * cos((1) * t);
        
%         nt = ones(m, 1) * ((20) * cos((1) * t) +...
%                             (10) * sin((5) * t) +...
%                              ((2) * cos((10) * t)) * exp(-0.2*t) );


    % ***********************
    %
    % RADP ON CART INVERTED PENDULUM SYSTEM -- SWEEP x_0
    %   
    
    case 'radp_cip_sweep_x0'
        
%         nt = ones(m, 1) * (20) * cos((1) * t);
        
        nt = ones(m, 1) * ((20) * sin((1) * t)...
                + (10) * sin((10) * t)...
                + (5) * cos(15 * t) );


    
    % ***********************
    %
    % COMPARISON D. VRABIE, F.L. LEWIS (2009) NONLINEAR "HARD" EXAMPLE
    % (SEC. 6.2)
    %   

    case 'comp_vrabie_lewis_2009_hard'
        
        nt = ones(m, 1) * (5) * sin((1) * t);
        
%         nt = ones(m, 1) * ((5) * sin((1) * t) +...
%                             (10) * sin((3) * t))    ;
                        
%         nt = 2 * ones(m, 1) * ((5) * sin((1) * t) +...
%                             (10) * sin((3) * t))+...
%                             ((2) * cos((4) * t)) * exp(-0.1*t);                        

    % ***********************
    %
    % COMPARISON OF CART INVERTED PENDULUM SYSTEM
    %   

    case 'comp_cip'
        
%         nt = ones(m, 1) * (50) * cos((1) * t);
        
        nt = ones(m, 1) * ( (50) * cos((1) * t)...
                            + (25) * sin(5 * t) );

%         nt = ones(m, 1) * ( (50) * cos((1) * t)...
%                             + (25) * sin(5 * t)...
%                             + (25) * sin(10 * t) );                        
                        
%         nt = ones(m, 1) * ((10) * sin((1) * t)...
%                 + (3) * sin((10) * t));

%         nt = 10 * ones(m, 1) * ((10) * cos((1) * t)...
%                 + (5) * sin((10) * t));

%         nt = 5 * ones(m, 1) * ((20) * sin((1) * t)...
%                 + (10) * sin((10) * t)...
%                 + (5) * cos(15 * t) );

%         nt = 10 * ones(m, 1) * ((10) * cos((1) * t)... 
%                 + (10) * sin((5) * t)...
%                 + (10) * sin((10) * t)...
%                 + (5) * cos(15 * t) );            

    % *********************************************************************
    % *********************************************************************
    %
    % THROW ERROR IF TAG DOES NOT COME UP A MATCH
    %   
    % *********************************************************************
    % *********************************************************************
    
    otherwise
        
        error('*** ERROR: NOISE TAG NOT RECOGNIZED ***');

end