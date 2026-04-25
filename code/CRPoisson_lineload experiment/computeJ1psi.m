function [j1] = computeJ1psi(s4e, pos4n, dirichlet_nodes, counts_per_node, sides)
% This function computes the coefficients of the J1 operator applied to the
% Crouzeix-Raviart basis functions for each given side.
% It averages the nodal values from neighboring triangles and returns a
% matrix of size size(counts_per_node,1) x size(sides,2). Boundary nodes 
% are explicitly set to zero.
%
% inputs: s4e             nr_elems x 3 matrix from computeS4e
%         pos4n           nr_nodes x 3nr_nodes matrix from computePos4n
%         n4sDb           nr_dirichletNodes x 2 matrix that gives the two nodes
%                         for each boundary edge
%         counts_per_node nr_nodes long vector with the count of occurences
%                         in n4e
%         sides           vector with integers, that give the side indeces
%                         (as in s4e) where the CR basis functions are 1
%
% output: j1              nr_nodes x length(sides) matirx with the
%                         coefficients for the J1 smoothing of a CR
%                         basis function in each column

nr_elems = size(s4e,1);
nr_sides = length(sides);

% Reorder matrix such that side j is opposite to node j 
% (local enumeration with j=1,2,3)
s4e = s4e(:, [2 3 1]);

% 1. Ensure sides is a row vector for implicit expansion (broadcasting)
sides = sides(:)'; 

% 2. Create sparse logical masks for matches across all elements and sides
% Each matrix M is of size (nr_elems x num_sides)
[tf1, loc1] = ismember(s4e(:,1), sides);
M1 = sparse(find(tf1), loc1(tf1), 1, nr_elems, nr_sides);

[tf2, loc2] = ismember(s4e(:,2), sides);
M2 = sparse(find(tf2), loc2(tf2), 1, nr_elems, nr_sides);

[tf3, loc3] = ismember(s4e(:,3), sides);
M3 = sparse(find(tf3), loc3(tf3), 1, nr_elems, nr_sides);

% Identify which elements contain the side at all (1 if yes, 0 if no)
HasMatch = M1 | M2 | M3;

% 3. Construct the Z components arithmetically 
% Base value is 1 (HasMatch), subtract 2 where the exact match occurs
Z1 = HasMatch - 2 * M1;
Z2 = HasMatch - 2 * M2;
Z3 = HasMatch - 2 * M3;

% 4. Assemble Z_long matching the behavior of Z(:) 
% Size becomes (3*nr_elems) x num_sides
Z_long_all = [Z1; Z2; Z3];

% 5. Sum Z contributions for each node via sparse matrix-vector multiplication
% Resulting size: (nr_nodes x num_sides)
sums_per_node = pos4n * Z_long_all;

% 6. Compute average (implicit expansion over columns)
j1 = sums_per_node ./ counts_per_node(:);

% 7. Enforce Dirichlet nodes to zero across all side columns at once
if ~isempty(dirichlet_nodes)
    j1(dirichlet_nodes, :) = 0;
end
end