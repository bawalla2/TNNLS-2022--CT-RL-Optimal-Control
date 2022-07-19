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

% Vector of x values to evaluate value function and policy at
x1vec = [-1 : 0.1 : 1];
x2vec = [-1 : 0.1 : 1];

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

%%
% *********************************************************************
%
% CONFIGURE BASES
%

% ***********************
%
% ORIGINAL BASES USED IN COMPARISON
%  

% Critic
critic_b1.tag = 'order_2_degree_4';   

c1 = [0.5; 0; 1; 0; 0; 0; 0; 1];    % Weight vals. after convergence

critic_b1 = config_basis(critic_b1, b_sett);


% Actor
actor_b1.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';

w1 = [-1; -2];      % Weight vals. after convergence

actor_b1 = config_basis(actor_b1, b_sett);


% % ***********************
% %
% % MODIFIED -- ALL EVEN MONOMIALS DEG. <= 4, CUSTOM ODD MONOMIALS DEG. <= 3
% %  
% 
% % Crititc
% 
% % All even monomials of total degree <= 4
% critic_b2.tag = 'monomial';
% critic_b2.type = 'even';   
% critic_b2.K = 4;
% 
% c2 = [ 1.3792
%     3.6754
%    -0.3214
%     0.7959
%     0.8907
%    -0.9242
%    -2.5847
%     0.1544 ];   % Weight vals. after convergence
% 
% 
% % Actor
% 
% % Manually-chosen odd monomials of total degree <= 3
% actor_b2.tag = 'monomial';
% actor_b2.type = 'custom';    
% actor_b2.dmat = [   
%                             1 0     % x_1
%                             0 1     % x_2
%                             1 2     % x_1 * x_2^2
%                             2 1     % x_1^2 * x_2
%                             3 0     % x_1^3
%                             0 3     % x_2^3
%                                 ];
% 
% w2 = [   -0.1871
%    -2.1764
%     3.1856
%    -0.8207
%     0.3786
%    -1.9610 ];   % Weight vals. after convergence                          


% ***********************
%
% MODIFIED -- DEBUG BASIS, OPTIMAL CRITIC BASIS
%  

% Crititc
critic_b2.tag = 'comp_vrabie_lewis_2009_hard_critic_debug';

c2 = [ -0.3235
   -0.3403
    0.2459
    0.1232
    0.5207
    0.0946
    0.9725
    1.6004 ];   % Weight vals. after convergence


% Actor
actor_b2.tag = 'comp_vrabie_lewis_2009_hard_actor_no_g';

w2 = [ -0.4457
   -3.4211 ];   % Weight vals. after convergence
    

% ***********************
%
% CONFIGURE MODIFIED BASES
% 

% Configure critic
critic_b2 = config_basis(critic_b2, b_sett);

% Configure actor
actor_b2 = config_basis(actor_b2, b_sett);                            
                            
%%
% *********************************************************************
%
% EVALUATE AND PLOT OPTIMAL VALUE FUNCTION ESTIMATES
%


nx1 = length(x1vec);
nx2 = length(x2vec);
[X,Y] = meshgrid(x1vec, x2vec);

V1x = zeros(nx1,nx2);
V2x = zeros(nx1,nx2);

for i = 1:nx1
   
    for j = 1:nx2
       
        x = [x1vec(i) ; x2vec(j)];
        
        [phi1x, ~] = eval_phi(x, critic_b1);
        [phi2x, ~] = eval_phi(x, critic_b2);
        
        V1x(i,j) = phi1x' * c1;
        V2x(i,j) = phi2x' * c2;
        
    end
    
end


% PLOT -- BOTH VALUE FUNCTIONS
figure(figcount)
h_fig = surf(X, Y, V1x, ...
    'FaceAlpha',0.5,'EdgeColor','red', 'EdgeAlpha',0.5);
hold on 
h_fig = surf(X, Y, V2x, 'FaceAlpha',0.5,'EdgeColor','none');
% set(h_fig, 'LineWidth', 2)
title('Cost Function')
grid on
legend('Optimal', 'Modified')
xlabel('x_1');
ylabel('x_2');
zlabel('V(x)');

% SAVE PLOT
if savefigs
    filename = ['V_x'];
    savepdf(figcount, relpath, filename); 
end

% Increment figure counter
figcount = figcount + 1; 

% PLOT -- VALUE FUNCTION DIFFERENCES
figure(figcount)
h_fig = surf(X, Y, V1x-V2x, ...
    'FaceAlpha',0.5,'EdgeColor','red', 'EdgeAlpha',0.5);
title('Cost Function Difference')
grid on
xlabel('x_1');
ylabel('x_2');
zlabel('V_{opt.}(x) - V_{mod.}(x)');

% SAVE PLOT
if savefigs
    filename = ['diff_V_x'];
    savepdf(figcount, relpath, filename); 
end

% Increment figure counter
figcount = figcount + 1; 


%%
% *********************************************************************
%
% EVALUATE AND PLOT OPTIMAL VALUE FUNCTION ESTIMATES
%

u1x = zeros(nx1,nx2);
u2x = zeros(nx1,nx2);

for i = 1:nx1
   
    for j = 1:nx2
       
        x = [x1vec(i) ; x2vec(j)];
        
        [phi1x, ~] = eval_phi(x, actor_b1);
        [phi2x, ~] = eval_phi(x, actor_b2);
        
        u1x(i,j) = phi1x' * w1;
        u2x(i,j) = phi2x' * w2;
        
    end
    
end


% PLOT
figure(figcount)
h_fig = surf(X, Y, u1x, ...
    'FaceAlpha',0.5,'EdgeColor','red', 'EdgeAlpha',0.5);
hold on 
h_fig = surf(X, Y, u2x, 'FaceAlpha',0.5,'EdgeColor','none');
% set(h_fig, 'LineWidth', 2)
title('Policy')
grid on
legend('Optimal', 'Modified')
xlabel('x_1');
ylabel('x_2');
zlabel('\mu(x)');

% SAVE PLOT
if savefigs
    filename = ['u_x'];
    savepdf(figcount, relpath, filename); 
end

% Increment figure counter
figcount = figcount + 1; 

% PLOT -- POLICY DIFFERENCES
figure(figcount)
h_fig = surf(X, Y, u1x-u2x, ...
    'FaceAlpha',0.5,'EdgeColor','red', 'EdgeAlpha',0.5);
title('Policy Difference')
grid on
xlabel('x_1');
ylabel('x_2');
zlabel('\mu_{opt.}(x) - \mu_{mod.}(x)');

% SAVE PLOT
if savefigs
    filename = ['diff_u_x'];
    savepdf(figcount, relpath, filename); 
end

% Increment figure counter
figcount = figcount + 1; 

