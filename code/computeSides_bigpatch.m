function sides_bigpatch = computeSides_bigpatch(n4e,e4s,e4n_logical,s4e,sides_line)
%% this function computes all the sides of the big patch of a given array
%% of sides. The big patch includes all the elements that share a node
%% with the triangles that have the sides from sides_line. 
% input     n4e             -   nodes for elements (nr_elems x 3)
%           e4s             -   elements for sides (nr_sides x 2)
%           e4n_logical     -   elements for nodes (nr_nodes x nr_elems)
%           s4e             -   sides for elements (nr_elements x 3)
%           sides_line      -   sides on the line (nr_sides_line x 1)
%
% output    sides_bigpatch  -   array of all unique sides in the big patch
elems_smallpatch = e4s(sides_line,:);

nodes_smallpatch = n4e(elems_smallpatch(:),:);
nodes_smallpatch = unique(nodes_smallpatch(:));

%elems_bigpatch = find(any(e4n(nodes_smallpatch), 1));
%sides_bigpatch = s4e(elems_bigpatch')
elems_bigpatch = any(e4n_logical(nodes_smallpatch,:),1);
sides_bigpatch = s4e(elems_bigpatch,:);
sides_bigpatch = unique(sides_bigpatch(:));

end