function sides_nodepatch = computeSides_nodepatch(nodes,pos4n,s4e)
%% This function computes an array of unique sides (idexed as in s4e)
%% that are contained the in the nodepatch around the nodes given in nodes.
%% For a given node, the nodepatch contains all the triangles that touch
%% this node. 
% input:  nodes                 vector containing the node indices for 
%                               which the nodepatch should be computes
%         pos4n                 (nr_nodes x 3*nr_elems) matrix from
%                               computePos4n
%         s4e                   (nr_elems x 3) matrix from computeS4e
%
% output: sides_nodepatch       vector containing the ideces of the sides
%                               (as indexed in s4e) that are part of the
%                               triangles in the nodepatch

% find all triangles in nodepatch
nr_elems = size(s4e,1);
[~,positions] = find(pos4n(nodes,:));
elements = mod(positions-1,nr_elems)+1;
elements = unique(elements);

% find all edges from those triangles
sides_nodepatch = s4e(elements,:);
sides_nodepatch = unique(sides_nodepatch(:));

end