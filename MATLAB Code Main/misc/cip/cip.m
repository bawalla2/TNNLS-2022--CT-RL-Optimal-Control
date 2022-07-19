% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CART INVERTED PENDULUM SYSTEM
%
% A.A. Rodriguez, Brent Wallace
%
% 12/9/2020
%
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

%%
% *************************************************************************
% *************************************************************************
%
% INITIALIZATION
% 
% *************************************************************************
% *************************************************************************

% ***********************
%
% FIGURES
%

% Close figures, clear command window.
close all
clc

% Figure saving controls.
savefigs = 0;           % Save figures control.
relpath = 'figures/';  	% Relative file path for saved figures.
figcount = 1;           % Initialize figure counter.

% Create save directory if figures are to be saved.
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
    relpath = [relpath, timestamp];   % Update relative path.
    mkdir(relpath);                   % Create directory for relative path.
end

% ***********************
%
% ADD SDToolbox
%
% addpath('../../../SDToolbox');


% ***********************
%
% PLOT SETTINGS
%

wmin = 1e-2;        % Minimum frequency for frequency response plots.
wmax = 1e2;         % Maximum frequency for frequency response plots.
numwpts = 1000;     % Number of frequency points for plots.

% Vector of frequency points.
wvec = logspace(log10(wmin),log10(wmax),numwpts);


% *************************************************************************
%
% SYSTEM PARAMETERS
% 

% % AAR VALUES
% mc = 0.455;             % Cart mass (kg).
% mp = 0.210;             % Pendulum mass (kg).
% l = 0.3048;             % Pendulum length (m).
% g = 9.8;                % Gravitational field constant (m/s^2).

% BAW VALUES
mc = 1;              % Cart mass (kg).
mp = 0.1;            % Pendulum mass (kg).
l = 0.5;             % Pendulum length (m).
g = 9.81;            % Gravitational field constant (m/s^2).

%
%****************************************
%
% Cart Inverted Pendulum Model
%
Ap = [  0       1       0               0
        0       0       -(mp*g/mc)        0
        0       0       0               1
        0       0       (g/l)*(1+(mp/mc)) 0  ];

Bp = [  0
        1/mc
        0
        -1/(mc*l) ];
    
Cp = [ 1   0   0   0];

Dp = 0;

P_ux = ss(Ap,Bp,Cp,Dp);


%%
% *************************************************************************
% *************************************************************************
%
% BEGIN MAIN
% 
% *************************************************************************
% *************************************************************************


%
%****************************************
%
% Check Controllability
%
conmat_cond_number = cond(ctrb(Ap,Bp))
conmat_rank = rank(ctrb(Ap,Bp))


% ***********************
%
% PLANT FREQUENCY RESPONSE
%

axes_vec_mag = [wmin,wmax,-40,40];      % Window to be used for magnitude
axes_vec_ph = [wmin,wmax,-200,0];      % Window to be used for phase

% See plotbode.m under the folder "SDToolbox" for documentation

sys_cell = {P_ux};
lgd_text = [];
wvec_cell = wvec;

plottypes = [1;1];  
ttl_cell = {'Plant P_{ux} Magnitude Response','Plant P_{ux} Phase Response'};
axes_cell = {axes_vec_mag, axes_vec_ph};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);

% SAVE FIGURE
filename = ['P_mag'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter

% SAVE FIGURE
filename = ['P_ph'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter

%%


%
%****************************************
%
% Desired Closed Loop Poles for Pole Placement Design
%
clpoles = [ -1 ;
            -2 ; 
            -0.5000+1.3229i;
            -0.5000-1.3229i ]   % NOTE: No eigenvalue should have multiplicity 
                                %       greater than number of inputs
                                
clpoles = [ -50 ;
            -75 ; 
           -0.5000+1.3229i;
            -0.5000-1.3229i ]   % NOTE: No eigenvalue should have multiplicity 
                                %       greater than number of inputs
%
%****************************************
%
% Pole Placement Design
%
[g1, prec, message] = place(Ap,Bp,clpoles) % PREC measures the number of accurate decimal digits 
                                        % in the actual closed-loop poles).  If some nonzero 
                                        % closed-loop pole is more than 10% off from the 
                                        % desired location, MESSAGE contains a warning message. 

final_clpoles1 = eig(Ap-Bp*g1)             % Final Closed Loop Poles



%   M. Wette 10-1-86
%   Revised 9-25-87 JNL
%   Revised 8-4-92 Wes Wang
%   Revised 10-5-93, 6-1-94 Andy Potvin
%   Revised 4-11-2001 John Glass, Pascal Gahinet
%
%   Ref:: Kautsky, Nichols, Van Dooren, "Robust Pole Assignment in Linear 
%         State Feedback," Intl. J. Control, 41(1985)5, pp 1129-1155


%  ACKER  Pole placement gain selection using Ackermann's formula.
% 
%   K = ACKER(A,B,P)  calculates the feedback gain matrix K such that
%    the single input system
%            .
%            x = Ax + Bu 
% 
%    with a feedback law of  u = -Kx  has closed loop poles at the 
%    values specified in vector P, i.e.,  P = eig(A-B*K).
% 
%    Note: This algorithm uses Ackermann's formula.  This method
%    is NOT numerically reliable and starts to break down rapidly
%    for problems of order greater than 10, or for weakly controllable
%    systems.  A warning message is printed if the nonzero closed-loop
%    poles are greater than 10% from the desired locations specified 
%    in P.


%
%****************************************
%
% Closed Loop Frequency Response
%
figure(1)
winit = -2;
wfin  = 2;
nfpts = 300;
w = logspace(winit,wfin,nfpts);
sys1 = ss(Ap-Bp*g1, Bp, g1, 0);
sv1 = sigma(sys1,w);
semilogx(w,20*log10(sv1))
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
title('Singular Values of G(sI - A + BG)^{-1} B')
%return
%
%****************************************
%
% Observability - ALL POSSIBILITIES (RANK = 4,3,2)
%
%
%**********************************************
% 
% RANK = 4 - POSITION INCLUDED
%
% 1
c = [ 1 0 0 0];                         % POSITION - OBSERVABLE - rank = 4, cond = 4.5231
[obsmat_cond_number_rank1 ] = [ cond(ctrb(Ap',c'))   rank(ctrb(Ap',c')) ]

% 12
c = [ 1 0 0 0
      0 1 0 0];                         % X, SPEED - NOT OBSERVABLE - rank = 4, cond = 212.5041
[obsmat_cond_number_rank12 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]

% 13
c = [ 1 0 0 0
      0 0 1 0];                         % X, THETA, - NOT OBSERVABLE - rank = 4, cond = 47.1888
[obsmat_cond_number_rank13 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]

% 14
c = [ 1 0 0 0
      0 0 0 1];                         % X, THETA DOT, - NOT OBSERVABLE - rank = 4, cond = 2205.8
[obsmat_cond_number_rank14 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]

% 123
c = [ 1 0 0 0
      0 1 0 0
      0 0 1 0];                         % X, THETA DOT, - NOT OBSERVABLE - rank = 4, cond = 217.633
[obsmat_cond_number_rank123 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]


% 124
c = [ 1 0 0 0
      0 1 0 0
      0 0 0 1];                         % X, THETA DOT, - NOT OBSERVABLE - rank = 4, cond = 2126.0
[obsmat_cond_number_rank124 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]

% 134
c = [ 1 0 0 0
      0 0 1 0
      0 0 0 1];                         % X, THETA DOT, - NOT OBSERVABLE - rank = 4, cond = 2206.3
[obsmat_cond_number_rank134 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]

% 1234
c = [ 1 0 0 0
      0 1 0 0
      0 0 1 0
      0 0 0 1];                         % X, THETA DOT, - NOT OBSERVABLE - rank = 4, cond = 2216.5
[obsmat_cond_number_rank1234 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]


%return
%
%**********************************************
% 
% RANK = 3 - POSITION NOT INCLUDED, SPEED INCLUDED
%
% 2
c = [ 0 1 0 0];                         % SPEED - NOT OBSERVABLE - rank = 3
[obsmat_cond_number_rank2 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]


% 23
c = [ 0 1 0 0
      0 0 1 0];                         % SPEED, THETA, - NOT OBSERVABLE - rank = 3
[obsmat_cond_number_rank23 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]

% 24
c = [ 0 1 0 0
      0 0 0 1];                         % SPEED, THETA DOT, - NOT OBSERVABLE - rank = 3
[obsmat_cond_number_rank24 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]

% 234
c = [ 0 1 0 0
      0 0 1 0
      0 0 0 1];                         % SPEED, THETA DOT, - NOT OBSERVABLE - rank = 3
[obsmat_cond_number_rank234 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]

%
%**********************************************
% 
% RANK = 2 - POSITION AND SPEED NOT INCLUDED
%
%
% 3
c = [ 0 0 1 0];                         % THETA - NOT OBSERVABLE - rank = 2
[obsmat_cond_number_rank3 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]


% 34
c = [ 0 0 1 0
      0 0 0 1];                         % THETA, THETA DOT - NOT OBSERVABLE - rank = 2
[obsmat_cond_number_rank34 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]


% 4
c = [ 0 0 0 1];                         % THETA DOT - NOT OBSERVABLE - rank = 2
[obsmat_cond_number_rank4 ] = [ cond(ctrb(Ap',c')) rank(ctrb(Ap',c'))]


%
%****************************************
%
% LQR DESIGNS
%
% LQR  Linear-quadratic regulator design for continuous-time systems.
%
%    [K,S,E] = LQR(A,B,Q,R,N)  calculates the optimal gain matrix K 
%    such that the state-feedback law  u = -Kx  minimizes the cost 
%    function
% 
%        J = Integral {x'Qx + u'Ru + 2*x'Nu} dt
%                                   .
%    subject to the state dynamics  x = Ax + Bu.
% 
%    The matrix N is set to zero when omitted.  Also returned are the
%    Riccati equation solution S and the closed-loop eigenvalues E:
%                        -1
%      SA + A'S - (SB+N)R  (B'S+N') + Q = 0 ,    E = EIG(A-B*K) .
% 
% Q = [   0.25    0       0       0
%         0       0       0       0
%         0       0       4       0
%         0       0       0       0 ]
    
% Q = [   .25     0       0       0
%         0       200     0       0
%         0       0       4       0
%         0       0       0       0 ]
   
% Q = [   0.25    0       0       0
%         0       0.2     0       0
%         0       0       4       0
%         0       0       0       0 ]

Q = 1 * eye(4)
    
%
%****************************************
%
% Closed Loop Frequency Response
%
figure(2)
r2 = 0.0001
[g2,k2,clpoles2] = lqr(Ap,Bp, Q,r2)
damp(Ap-Bp*g2)
w = logspace(winit,wfin,nfpts);
sys2 = ss(Ap-Bp*g2, Bp, g2, 0);
sv2 = sigma(sys2,w);
semilogx(w,20*log10(sv2))
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
title('Singular Values of G(sI - A + BG)^{-1} B')

figure(3)
r3 = 0.001
[g3,k3,clpoles3] = lqr(Ap,Bp, Q,r3)
damp(Ap-Bp*g3)
w = logspace(winit,wfin,nfpts);
sys3 = ss(Ap-Bp*g3, Bp, g3, 0);
sv3 = sigma(sys3,w);
semilogx(w,20*log10(sv3))
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
title('Singular Values of G(sI - A + BG)^{-1} B')

figure(4)
r4 = 0.01
[g4,k4,clpoles4] = lqr(Ap,Bp, Q,r4)
damp(Ap-Bp*g4)
w = logspace(winit,wfin,nfpts);
sys4 = ss(Ap-Bp*g4, Bp, g4, 0);
sv4 = sigma(sys4,w);
semilogx(w,20*log10(sv4))
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
title('Singular Values of G(sI - A + BG)^{-1} B')

figure(5)
r5 = 0.1
[g5,k5,clpoles5] = lqr(Ap,Bp, Q,r5)
damp(Ap-Bp*g5)
w = logspace(winit,wfin,nfpts);
sys5 = ss(Ap-Bp*g5, Bp, g5, 0);
sv5 = sigma(sys5,w);
semilogx(w,20*log10(sv5))
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
title('Singular Values of G(sI - A + BG)^{-1} B')

figure(6)
r6 = 100
[g6,k6,clpoles6] = lqr(Ap,Bp, Q,r6)
damp(Ap-Bp*g6)
w = logspace(winit,wfin,nfpts);
sys6 = ss(Ap-Bp*g6, Bp, g6, 0);
sv6 = sigma(sys6,w);
semilogx(w,20*log10(sv6))
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
title('Singular Values of G(sI - A + BG)^{-1} B')


%%
%****************************************
%
% LQR DESIGNS -- STABILIZING (A, -1/2 B R^{-1} B^T)
%

R = 1;              % R as defined in the original problem


R_p = eye(4);       % "R" for the LQR problem associated with (A, -1/2 B R^{-1} B^T)
Q_p = eye(4);       % "Q" for the LQR problem associated with (A, -1/2 B R^{-1} B^T)

[K,G,clps] = lqr(Ap,-0.5*Bp*inv(R)*Bp', Q_p,R_p)

damp(Ap-0.5*Bp*inv(R)*Bp'*K)


%%
%****************************************
%
% BRUTE FORCE SEARCH FOR STABILIZING (A, -1/2 B R^{-1} B^T)
%

R = 1;              % R as defined in the original problem

c2vec = (-1000:1:1000)';
nc2 = size(c2vec,1);
c4vec = (-1000:1:1000)';
nc4 = size(c4vec,1);

for i = 1:nc2
    c2 = c2vec(i);
    for j = 1:nc4
        c4 = c4vec(j);
        Cij = diag([0,c2,0,c4]);
        Aij = Ap -0.5*Bp*inv(R)*Bp' * Cij;
        eigs = eig(Aij);
%         eigs
        maxreeig = max(real(eigs));
        if maxreeig <= 0
           disp(['Hurwitz for (c2, c4) = ('  num2str(c2) ', ' num2str(c4) ')'])
        end
    end    
end
