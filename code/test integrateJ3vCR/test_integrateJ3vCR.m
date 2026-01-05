%% This module tests the integrateJ3vCR function.

% import AFEM library
addpath(genpath('C:\Users\natha\Documents\MATLAB\afem-base-master'))

%% define test_case and error tolerance
% case 1 is an all zeros input
% case 2 is an all ones input
% case 3 is a conform CR function input
% case 4 is a single CR-basis function input
test_case = 1;
epsilon = 0.0001;

% load triangulation of a square with 32 triangles
[c4n, n4e, n4sDb] = loadGeometry('Square', 4);

% compute helpful matrices and constants
n4s = computeN4s(n4e);
s4n = computeS4n(n4e);
s4e = computeS4e(n4e);
area4e = computeArea4e(c4n,n4e);

nr_sides = size(n4s,1);
nr_vertices = size(c4n,1);

% get Direchlet boundary sides
DbSides = zeros(1,size(n4sDb,1));
for i = 1:size(n4sDb,1)
    DbSides(i) = s4n(n4sDb(i,1),n4sDb(i,2));
end


% define the test-input
switch test_case
    case 1
        vCR = zeros(nr_sides,1);
    case 2
        vCR = ones(nr_sides,1);
        vCR(DbSides) = 0;
    case 3
        x = (1:25)';
        vCR = computeCRfromP1(n4s,x);
        vCR(DbSides) = 0;
    case 4
        vCR = zeros(nr_sides,1);
        vCR(2) = 1;
        
end


% get coefficients for J3vCR
averaging_coefficients = computeJ1vCR(c4n,n4e,n4sDb,vCR);
bubble_coefficients = computeJ2vCR(n4s,vCR,averaging_coefficients);
volume_coefficients = computeJ3vCR(c4n,n4e,vCR,averaging_coefficients,bubble_coefficients);

% define the integrand J3vCR(lambda1, lambda2, lambda3)
integrand = @(n4p,Gpts4p,Gpts4ref) ...
            J3vCR(n4p,Gpts4p,Gpts4ref,n4e,s4e,averaging_coefficients,bubble_coefficients,volume_coefficients,area4e);

% compute exact and approximated integral
exact_integral = integrateJ3vCR(averaging_coefficients, bubble_coefficients, volume_coefficients, c4n, n4e);
approx_integral = sum(integrate(c4n,n4e,integrand,4,area4e));

% check if the results match
check = ((exact_integral <= approx_integral + epsilon) | (exact_integral >= approx_integral - epsilon))
