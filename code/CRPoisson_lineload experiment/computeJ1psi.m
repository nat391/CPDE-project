function [j1] = computeJ1psi(s4e,pos4n,n4sDb,counts_per_node,side)
% this function computes the coefficients of the J1 operator applied to a
% Crouzeix-Raviart function.
% It averages the nodal values from neighboring triangles and returns a
% vector of length size(c4n,1. Boundary nodes are explicitly set to zero.
% inputs: s4e               nr_elems x 3 matrix from computeS4e
%         pos4n             nr_nodes x 3nr_nodes matrix from computePos4n
%         n4sDb             |nr_dirichletNodes| x 2 matrix that gives the
%                           two nodes for each boundary edge
%         counts_per_node   nr_nodes long vector with the count of 
%                           occurences in n4e
%         side              integer which is the side index (as in s4e)
%                           where the CR basis function is 1
%
% output: j1                nr_nodes long vector with the coefficients
%                           of the J1 operator

% reorder matrix such that side j is opposite to node j 
% (local enumeration with j=1,2,3)
s4e = s4e(:,[2 3 1]);

% compute the node values Z for each Crouzeix-Raviart patch by extrapolating
% the given edge-values using: psi_j = 1 - 2lambda_j (see GNUMO p.157)
[i,j] = find(s4e == side);

nr_elems = size(s4e,1);
Z = sparse(nr_elems, 3);

for k = 1:numel(i)
    Z(i(k),:) = [1 1 1];
    Z(i(k),j(k)) = -1;
end

% Flatten Z so positions align
Z_long = Z(:);      % same ordering

% Sum Z contributions for each node (vectorized via sparse mat-vec)
% we can do: sums(node) = sum of Z_long at positions where node occurs
sums_per_node = pos4n * Z_long;        % (nr_nodes x 1)

% Compute average
j1 = sums_per_node ./ counts_per_node;


% Enforce Dirichlet nodes to zero (if required)
if ~isempty(n4sDb)
    j1(unique(n4sDb(:))) = 0;
end

