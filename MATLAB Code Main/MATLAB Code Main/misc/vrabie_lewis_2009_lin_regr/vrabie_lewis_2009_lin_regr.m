%%
% *********************************************************************
%
% CONFIG
%

% Close all
close all
clear
clc

% ADD PATHS
addpath('..\..\eval_functs');
addpath('..\..\config');

% SYS AND ALG
sys.tag = 'vrabie_lewis_2009_hard';
alg = 'debug';

[sys, sys_plot_settings] = config_sys(sys);

% Settings for basis initialization
b_sett.sys = sys;
b_sett.alg = alg;

% Grid of x values to perform regression over
x1vec = [-1 : 0.1 : 1]';
x2vec = [-1 : 0.1 : 1]';
% x1vec = [-1 : 0.25 : 1]';
% x2vec = [-1 : 0.25 : 1]';
% x1vec = [-1 : 0.5 : 1]';
% x2vec = [-1 : 0.5 : 1]';

% x1vec = [-1 : 0.1 : 1]';
% x2vec = [0 : 0.01 : 1]';

% Grid of x values to plot regression function over
x1pvec = [-1 : 0.25 : 1]';
x2pvec = [-1 : 0.25 : 1]';

% Save figures or not
savefigs = 0;

% Initialize figure counter
figcount = 1;

% Relative path to figures folder
relpath = '..\..\00 figures';

% Create save directory if figures are to be saved
if savefigs
    time = fix(clock);
    timestamp = '';
    for i = 1:length(time)
       timestamp = [timestamp , num2str(time(i))];
       if i < length(time)
           timestamp = [timestamp , '_'];
       end
    end
    timestamp = [timestamp, '\'];
    relpath = [relpath, timestamp];   % Update relative path
    mkdir(relpath);                   % Create directory for relative path
end


% ***********************
%
% SURFACE PLOT SETTINGS -- WEIGHT VALUES
%        

% Surface plot face color
facecolor = 'blue';

% Surface plot face transparency
facealpha = 0.25;

% Surface plot edge transparency
edgealpha = 0.75;

% ***********************
%
% SURFACE PLOT SETTINGS -- OPTIMAL WEIGHT VALUE PLANES
%        

% Surface plot face color
facecolor_opt = 'red';

% Surface plot face transparency
facealpha_opt = 0.1;

% Surface plot edge transparency
edgealpha_opt = 0.5;       



%%
% *********************************************************************
%
% CONFIGURE BASES
%


% ***********************
%
% OPTIMAL ACTOR BASIS
%   

% Optimal actor basis
b_mu_opt.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';
b_mu_opt = config_basis(b_mu_opt, b_sett);

% Optimal actor basis weights
w_star = [-1 ; -2];


% ***********************
%
% ACTOR BASIS TO PERFORM REGRESSION ON
%   

% All odd monomials of total degree <= 2
b_mu_hat.tag = 'monomial';
b_mu_hat.type = 'even';   
b_mu_hat.K = 4;

b_mu_hat = config_basis(b_mu_hat, b_sett);
                   
                            
%%
% *********************************************************************
%
% EVALUATE FUNCTIONS, PERFORM REGRESSION
%

% ***********************
%
% EVALUATE FUNCTIONS
%  

nx1 = length(x1vec);
nx2 = length(x2vec);
[X1,X2] = meshgrid(x1vec, x2vec);

% Data storage
mu_star_x = zeros(nx1*nx2, 1);
psi_x_mat = zeros(nx1*nx2, b_mu_hat.N);

% Initialize counter
count = 1;

% Evaluate functions
for j = 1:nx2
   
    for i = 1:nx1
       
        % Extract x
        x = [x1vec(i) ; x2vec(j)];
        
        % Evaluate bases
        psix_mustar = eval_phi(x, b_mu_opt);
        psix_muhat = eval_phi(x, b_mu_hat);
        
        % Evaluate optimal policy \mu^*(x)
        mustar_x = w_star' * psix_mustar;
        
        % Store \mu^*(x), \psi(x)
        mu_star_x(count) = mustar_x;
        psi_x_mat(count,:) = psix_muhat';
        
        % Increment counter
        count = count + 1;
        
    end
    
end


% ***********************
%
% PERFORM REGRESSION
%  

% Regression
w_lsq = psi_x_mat \ mu_star_x;

w_lsqp = lsqr(psi_x_mat, mu_star_x);

% Approximate optimal policy
mu_lsq_x = psi_x_mat * w_lsq;


%%
% *********************************************************************
%
% PLOT RESULTS OF REGRESSION
%

% ***********************
%
% EVALUATE FUNCTIONS
%  

nx1p = length(x1pvec);
nx2p = length(x2pvec);
[X1p,X2p] = meshgrid(x1pvec, x2pvec);

% Data storage
mu_star_xp = zeros(nx1p*nx2p, 1);
mu_lsq_xp = zeros(nx1p*nx2p, 1);

% Initialize counter
count = 1;

% Evaluate functions
for j = 1:nx2p
   
    for i = 1:nx1p
       
        % Extract x
        x = [x1pvec(i) ; x2pvec(j)];
        
        % Evaluate bases
        psix_mustar = eval_phi(x, b_mu_opt);
        psix_muhat = eval_phi(x, b_mu_hat);
        
        % Evaluate optimal policy \mu^*(x)
        mustar_x = w_star' * psix_mustar;
        
        % Evaluate least-squares optimal policy
        mulsq_x = w_lsq' * psix_muhat;
        
        % Store \mu^*(x), \psi(x)
        mu_star_xp(count) = mustar_x;
        mu_lsq_xp(count) = mulsq_x;
        
        % Increment counter
        count = count + 1;
        
    end
    
end

% Reshape data from nx1*nx2 x 1 to nx1 x nx2
mu_star_xp = reshape(mu_star_xp, [nx1p,nx2p]);
mu_lsq_xp = reshape(mu_lsq_xp, [nx1p,nx2p]);

% PLOT -- OPTIMAL POLICY AND ESTIMATE OPTIMAL POLICY
figure(figcount)
h_fig = surf(X1p, X2p, mu_star_xp');
set(h_fig, 'FaceColor', facecolor_opt);
set(h_fig, 'FaceAlpha', facealpha_opt);
set(h_fig, 'EdgeAlpha', edgealpha);
hold on
h_fig = surf(X1p, X2p, mu_lsq_xp');
set(h_fig, 'FaceColor', facecolor);
set(h_fig, 'FaceAlpha', facealpha);
set(h_fig, 'EdgeAlpha', edgealpha);            
ttl = ['Optimal Policy $\mu^*$ and Least-Squares Regr. Policy '...
        '$\hat{\mu}_{i^*}$'];
title(ttl, 'interpreter', 'latex')
grid on
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
zlabel('$\mu(x)$', 'interpreter', 'latex');
lgd = legend({'$\mu^*$', '$\hat{\mu}$'},...
        'interpreter', 'latex', ...
        'location', 'northeast');
%             set(lgd, 'location', 'best')

% SAVE PLOT
if savefigs
    filename = ['mu_x'];
    savepdf(figcount, relpath, filename); 
end

% Increment figure counter
figcount = figcount + 1;   

