function [integral_values] = integrateJ3psi_lineload(c4n,n4s_line,sides_line,j1,j2,lineload,degree)
%% This function computes the exact integral of J3vCR along a line
%   input: 
%   c4n - coordinates for node
%   n4s_line   - nodes for sides along the line
%   sides_line - [nr_sides_line x 1] sides that lie on the line (indexing according to n4s)
%   j1         - (nr_nodes,k) sparse matrix that contains in each
%                column the coefficients of J1 applied to a CR function. 
%   j2         - (nr_sides,k) spars matrix that contains in each 
%                column the coefficients of J2 applied to a CR function.
%   lineload   - function handle, function by which to multiply J3vCR in 
%                the line integral
%   output: 
%   integral_values - (k,1) matrix containing the values of the
%                     integrals along the line for each given functions.

%% comments
%  1. the volume coefficients are not needed because they are zero along
%  the edges
%  2. we assume that the mesh "matches" the line

nr_sides_line = length(sides_line);

% define the integrand J3vCR(lambda1, lambda2, lambda3)
% by multiplication of sparse (s,s) and sparse (s,k) matrix
integrand = @(n4s_line, Gpts4p, Gpts4ref) ...
    spdiags(lineload(Gpts4p), 0, nr_sides_line, nr_sides_line) * ...
    J3vCR_line(n4s_line, Gpts4p, Gpts4ref, sides_line, j1, j2);

% get the integral value for each side and each function
val_matrix = integrate_sparse(c4n, n4s_line, integrand, degree);

% sum over all sides for each function and transpose
integral_values = sum(val_matrix, 1)';