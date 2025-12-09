%test

% import AFEM library
addpath(genpath('C:\Users\natha\Documents\MATLAB\afem-base-master'))

[c4n, n4e, n4sDb] = loadGeometry('Square', 2);
n4s = computeN4s(n4e);

vCR = zeros(size(n4s,1),1);
vCR(2) = 1;

averaging_coefficients = computeJ1vCR(c4n,n4e,n4sDb,vCR);

bubble_coefficients = computeJ2vCR(n4s,vCR,averaging_coefficients)
