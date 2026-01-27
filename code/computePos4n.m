function [pos4n,counts_per_node] = computePos4n(n4e,nr_nodes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nr_elements = size(n4e,1);
nr_positions = 3*nr_elements;

node_idx = n4e(:);    % (nr_positions x 1) node index of each position
pos_idx = (1:nr_positions)';  % column indices (positions)

% Create a sparse matrix with 1 where node occurs at a given position
pos4n = sparse(node_idx, pos_idx, 1, nr_nodes, nr_positions);  % (nr_nodes x nr_positions)

% Count occurrences per node
counts_per_node = sum(pos4n, 2); % (nr_nodes x 1)
end
