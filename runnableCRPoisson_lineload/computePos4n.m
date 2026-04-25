function [pos4n,counts_per_node] = computePos4n(n4e,nr_nodes)
% This function computes the (nr_nodes x 3nr_nodes) boolean matrix pos4n
% and the nr_nodes long vector counts_per_node. pos4n(i,j) = 1 if n4e(j)=i 
% and 0 otherwise. counts_per_node(i) gives the count of occurences of the
% node i in n4e. 
% input:    n4e             nr_elems x 3 matrix from loadGeometry
%           nr_nodes        integer = size(c4n,1)
%
% output:   pos4n           (nr_nodes x 3nr_nodes) matrix with pos4n(i,j) = 1
%                           if n4e(j)=i and 0 otherwise
%           counts_per_node nr_nodes long vector where counts_per_node(i) 
%                           gives the count of occurences of the node i in
%                           n4e.

nr_elements = size(n4e,1);
nr_positions = 3*nr_elements;

node_idx = n4e(:);    % (nr_positions x 1) node index of each position
pos_idx = (1:nr_positions)';  % column indices (positions)

% Create a sparse matrix with 1 where node occurs at a given position
pos4n = sparse(node_idx, pos_idx, 1, nr_nodes, nr_positions);  % (nr_nodes x nr_positions)

% Count occurrences per node
counts_per_node = full(sum(pos4n, 2)); % (nr_nodes x 1)
end
