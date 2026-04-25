function [j2] = computeJ2psi(n4s,sides,j1)
% this function computes the coefficients j2 of the J2 operator applied 
% to CR basis functions. 
% input:  n4s             nr_sides x 2 matrix from computeN4s
%         sides           vector containing the side indeces where CR basis 
%                         functions are 1 
%         j1              nr_nodes x length(sides) matrix from computeJ1psi
%
% output: j2              nr_sides x length(sides) matrix with the
%                         coefficients of the J2 smoothing applied for each 
%                         CR basis function defined by the indices in sides.

nr_sides = size(n4s, 1);
nr_selected_sides = length(sides);

% Ensure sides is a column vector for sparse matrix indexing
sides = sides(:);

% 1. Extract the node arrays for the two ends of each side
node1 = n4s(:, 1);
node2 = n4s(:, 2);

% 2. Compute integral mean of (vCR - J1psi) over edge F using midpoint rule
% This extracts (nr_sides x num_sides) matrices and computes the average
j2 = -(j1(node1, :) + j1(node2, :)) / 2;

% 3. Vectorize the +1 correction at j2(side(k), k)
% Construct a sparse matrix with 1s at (sides(k), k) for all k
correction = sparse(sides, (1:nr_selected_sides)', 1, nr_sides, nr_selected_sides);

% Apply the correction globally
j2 = j2 + correction;

end
