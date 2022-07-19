function c0 = config_basis_ICs(dmat, dmat_ic, c0_ic)
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
% This program sets the ICs for a monomial basis which was initialized
% automatically (i.e., of type 'all', 'even', or 'odd', cf.
% eval_phi.m).
%
% *************************************************************************
%
% CALL SYNTAX
%
% *************************************************************************
%
% c0 = config_basis_ICs(dmat, dmat_ic, c0_ic)
%
% *************************************************************************
%
% INPUTS
%
% *************************************************************************
% 
% dmat          (N x n Matrix) Matrix of monomial coefficients. See
%               eval_phi.m for full description.
% dmat_ic       ('numIC' x n Matrix) A matrix of 'numIC' <= N rows (where
%               'numIC' denotes the number of network weights to
%               initialize) and n columns, whose i-th row contains the
%               monomial coefficients corresponding to the desired monomial
%               to have the initial weight c0_ic(i). 'dmat' will be
%               searched row-wise for each row of 'dmat_ic' to find find
%               the position of the respective monomial in the basis.
% c0_ic         ('numIC'-dimensional vector) A vector consisting of <= N
%               elements, whose i-th entry corresponds to the initial
%               condition of the respective monomial in the i-th row of
%               'dmat_ic'.
%
% *************************************************************************
%
% OUTPUTS
%
% *************************************************************************
%
% c0            (N-dim. Vector) Vector of weight ICs matching 'c0_ic' in
%               the respective positions, and is zero in its other entries.
%
% *************************************************************************
% *************************************************************************
% *************************************************************************



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

% Number of basis functions
N = size(dmat,1);
numIC = size(dmat_ic,1);    % Number of entries of c0 to set

% Initialize IC vector
c0 = zeros(N,1);

for i = 1:numIC
  
    % Find the row of dmat with the monomial corresponding to the i-th row
    % of dmat_ic
    ind = find(ismember(dmat, dmat_ic(i,:), 'rows'));
    
    % Set c0(ind) to the corresponding value in c0_ic
    c0(ind) = c0_ic(i);
    
end



