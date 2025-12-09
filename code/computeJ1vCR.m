function [averaging_coefficients] = computeJ1vCR(c4n,n4e,n4sDb,vCR)
% this function computes the J1 smoothing of a Crouzeix-Raviart function.
% It averages the nodal values from neighboring triangles and returns a
% vector of length size(c4n,1)

% compute necessary helpers
nr_nodes = size(c4n,1); % remark: this is the only line where c4n is used, maybe pass the number of elements as an input instead?
nr_elements = size(n4e,1);
nr_positions = 3*nr_elements;
s4e = computeS4e(n4e);
% reorder matrix such that side j is opposite to node j 
% (local enumeration with j=1,2,3)
s4e = s4e(:,[2 3 1]);

% compute the node values Z for each Crouzeix-Raviart patch by extrapolating
% the given edge-values using: psi_j = 1 - 2lambda_j (see GNUMO p.157)
Z = sum(vCR(s4e),2)*ones(1,3)- 2*vCR(s4e);

%{
% compute the linear position of every node in the matrix n4e
nr_positions = 3*nr_elements;
p4n = zeros(nr_nodes, nr_positions);
for position = 1:nr_positions
    p4n(n4e(position), position) = 1;
end
p4n = logical(p4n);

% initialize solution vector averaging_coefficients
averaging_coefficients = zeros(nr_nodes,1);

% compute degrees of freedom (innner nodes)
dof=setdiff(1:nr_nodes,n4sDb(:));



% compute the average of the extrapolated node values Z for every inner node
for node = dof

    % find all indices at which the node occurs in the connectivity matrix
    %indices = 1:nr_positions(p4n(node,:));
    count = sum(p4n(node,:));
    
    % compute the average of the values at that node
    averaging_coefficients(node) = sum(Z(p4n(node,:)))/count;
end
end
%}
% Flatten connectivity and Z so positions align
node_idx = n4e(:);    % (nr_positions x 1) node index of each position (col-major)
Z_long   = Z(:);      % same ordering

% Build sparse position incidence: rows = node, cols = position
pos_idx = (1:nr_positions)';  % column indices (positions)
% Create a sparse matrix with 1 where node occurs at a given position
p4n = sparse(node_idx, pos_idx, 1, nr_nodes, nr_positions);  % (nr_nodes x nr_positions)
% convert to logical when needed:
p4n_logical = logical(p4n);

% Sum Z contributions for each node (vectorized via sparse mat-vec)
% we can do: sums(node) = sum of Z_long at positions where node occurs
sums_per_node = p4n * Z_long;        % (nr_nodes x 1)

% Count occurrences per node
counts_per_node = sum(p4n, 2);      % (nr_nodes x 1) (sparse sum is efficient)

% Avoid divide-by-zero: set zeros to 1 to keep 0/1 = 0, but better to mark them separately
zero_count_mask = (counts_per_node == 0);
counts_per_node(zero_count_mask) = 1;  % temporary safe divisor

% Compute average
averaging_coefficients = sums_per_node ./ counts_per_node;

% Set nodes with zero count to zero (if any)
averaging_coefficients(zero_count_mask) = 0;

% Enforce Dirichlet nodes to zero (if required)
if ~isempty(n4sDb)
    averaging_coefficients(unique(n4sDb(:))) = 0;
end

