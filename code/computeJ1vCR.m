function [averaging_coefficients] = computeJ1vCR(s4e,pos4n,n4sDb,counts_per_node,vCR)
% this function computes the J1 smoothing of a Crouzeix-Raviart function.
% It averages the nodal values from neighboring triangles and returns a
% vector of length size(c4n,1)

% reorder matrix such that side j is opposite to node j 
% (local enumeration with j=1,2,3)
s4e = s4e(:,[2 3 1]);

% compute the node values Z for each Crouzeix-Raviart patch by extrapolating
% the given edge-values using: psi_j = 1 - 2lambda_j (see GNUMO p.157)
Z = sum(vCR(s4e),2)*ones(1,3)- 2*vCR(s4e);

% Flatten Z so positions align
Z_long   = Z(:);      % same ordering

% Sum Z contributions for each node (vectorized via sparse mat-vec)
% we can do: sums(node) = sum of Z_long at positions where node occurs
sums_per_node = pos4n * Z_long;        % (nr_nodes x 1)

% Compute average
averaging_coefficients = sums_per_node ./ counts_per_node;


% Enforce Dirichlet nodes to zero (if required)
if ~isempty(n4sDb)
    averaging_coefficients(unique(n4sDb(:))) = 0;
end

