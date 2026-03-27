%test J1vCR

[c4n, n4e, n4sDb] = loadGeometry('square', 2);

s4e = computeS4e(n4e);
[pos4n, counts_per_node] = computePos4n(n4e,size(c4n,1));
n4s = computeN4s(n4e);
nr_sides = size(n4s,1);
vCR = zeros(nr_sides,1);
vCR(2) = 1;
%inner_sides = setdiff(1:nr_sides, n4sDb);
%vCR(inner_sides) = 1;

averaging_coefficients = computeJ1vCR(s4e,pos4n,n4sDb,counts_per_node,vCR);