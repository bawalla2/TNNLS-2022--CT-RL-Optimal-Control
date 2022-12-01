%% CONFIG

% ADD PATHS
addpath('..\..\eval_functs');
addpath('..\..\config');

% SYS AND ALG
sys.tag = 'vrabie_lewis_2009_hard';
alg = 'debug';

% SETTINGS FOR BASIS INITIALIZATION
b_sett.sys = sys;
b_sett.alg = alg;


[sys, sys_plot_settings] = config_sys(sys);

%% CONFIGURE BASES

b1.tag = 'order_2_degree_4';
        
b2.tag = 'monomial';
b2.type = 'custom';   
b2.dmat = [     
                            2 0
                            1 1
                            0 2
                            4 0
                            3 1
                            2 2
                            1 3
                            0 4
                                    ];

                                       
b_sett.sys = sys;
b_sett.alg = alg;

b1 = config_basis(b1, b_sett);
b2 = config_basis(b2, b_sett);

%% EVALUATE AT SINGLE POINT


% x = [14.3 ; -23.4];
x = [4975 ; -239.54085];
% x = [-39.304 ; 671.096];

[phix1, dphix1] = eval_phi(x, b1)
[phix2, dphix2] = eval_phi(x, b2)

norm(phix1 - phix2)
norm(dphix1(:) - dphix2(:))


%% EVALUATE AT A BUNCH OF RANDOM POINTS

nump = 5000;
mag = 1;
% mag = 100;
randx = mag * (rand(2,nump) - 0.5);
cumerr_phix = 0;
cumerr_dphix = 0;

for i = 1:nump
    x = randx(:,i);
    [phix1, dphix1] = eval_phi(x, b1);
    [phix2, dphix2] = eval_phi(x, b2);
    cumerr_phix = cumerr_phix + norm(phix1 - phix2);
    cumerr_dphix = cumerr_dphix + norm(dphix1(:) - dphix2(:)); 
end

cumerr_phix
cumerr_dphix
