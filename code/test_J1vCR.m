%test J1vCR
% import AFEM library
addpath(genpath('C:\Users\natha\Documents\MATLAB\afem-base-master'))

[c4n, n4e, n4sDb] = loadGeometry('square', 5);

n4s = computeN4s(n4e);
nr_sides = size(n4s,1);
vCR = zeros(nr_sides,1);
inner_sides = setdiff(1:nr_sides, n4sDb);
vCR(inner_sides) = 1;

averaging_coefficients = computeJ1vCR(c4n,n4e,n4sDb,vCR);