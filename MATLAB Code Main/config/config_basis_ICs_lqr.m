function w0 = config_basis_ICs_lqr(sys, dmat, Q, R)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% INITIALIZE BASIS PARAMETERS
%
% Brent Wallace  
%
% 2022-02-08
%
% This program finds the initial weights w \in R^{N_1} such that the critic
% network activation functions \Phi(x) satisfy:
%
%           \nabla \Phi^T(x) w = 2 * P x
%
% where P = P^T > 0 is the unique positive definite solution to the Riccati
% equation
%
%       0 = A^T P + P A - P B R^{-1} B^T P + Q
%
% and (A, B) denote the linearization of the nonlinear system (f, g). It
% can be shown that such weights w \in R^{N_1}, when implemented in the
% actor structure:
%
%       \hat{\mu}(x) = - 1/2 R^{-1} g^T(x) \nabla \Phi^T(x) w
%
% render the closed-loop nonlinear system locally asymptotically stable.
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% w0 = config_basis_ICs_lqr(sys, dmat, Q, R)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
% 
% sys           (Struct) Contains system parameters.
% dmat          (N_1 x n matrix) Matrix of monomial degree coefficients.
%               See eval_phi.m for full description. NOTE: This matrix is
%               always initialized programmatically by config_basis.m.
% Q             (n x n matrix) State penalty matrix Q = Q^T >= 0.
% R             (m x m matrix) Control penalty matrix R = R^T > 0.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% w0            (N_1-dim. Vector) Vector of weights satisfying the
%               properties described above.
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

% Extract system parameters
A = sys.A;          % System linearization "A" matrix
B = sys.B;          % System linearization "B" matrix

% System dimensions
n = size(A, 1);     % Number of states
m = size(B, 2);     % Number of inputs


% Number of basis functions
N = size(dmat,1);




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
% SOLVE LQR PROBLEM FOR LINEARIZATION
%
[K,P] = lqr(A,B,Q,R);

% Number of weight coefficients to set
numcoeff = n * (n + 1) / 2;

% Data storage
dmat_ic = zeros(numcoeff, n);
c0_ic = zeros(numcoeff, 1);

% ***********************
%       
% DETERMINE DESIRED WEIGHT COEFFICIENTS AND RESPECTIVE MONOMIALS
%

% Begin an accumulative counter
count = 1;

for i = 1:n
    
    for j = i:n
        
        % Determine the degree coefficients for the activation function x_i
        % * x_j
        if i == j
            dmat_ic(count, i) = 2;
        else
            dmat_ic(count, [i j]) = 1;
        end
        
        % Set the corresponding coefficient
        if i == j
            c0_ic(count) = P(i,j);
        else
            c0_ic(count) = 2 * P(i,j);
        end
       
        % Increment the counter
        count = count + 1;
        
    end
    
end
    

% ***********************
%       
% INITIALIZE THE WEIGHTS ACCORDING TO THE BASIS
%

w0 = config_basis_ICs(dmat, dmat_ic, c0_ic);


