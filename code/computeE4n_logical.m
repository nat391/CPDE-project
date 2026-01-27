function e4n_logical = computeE4n_logical(n4e, nr_nodes)
%% compute the (nr_nodes x nr_elements) matrix where 
%% e4n_logical(nodenr,elementnr) = 1 iff that node is part of that element
nr_elements = size(n4e,1);
node_idx = n4e(:);
element_idx = repelem((1:nr_elements)',3);
e4n_logical = sparse(node_idx,element_idx,1,nr_nodes,nr_elements);
end