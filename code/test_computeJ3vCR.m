% This module tests the computeJ3CR function. The calculations for 
% test_case 4 are outlined in "testingJ3_example.pdf"

%define test_case
test_case = 4;

% import AFEM library
addpath(genpath('C:\Users\natha\Documents\MATLAB\afem-base-master'))

% load triangulation of a square with 25 triangles
[c4n, n4e, n4sDb] = loadGeometry('Square', 2);

%compute the (2 x nr_sides) nodes for sides matrix
n4s = computeN4s(n4e);
nrSides = size(n4s,1);
nrNodes = size(c4n,1);

switch test_case
    case 1
        vCR = zeros(size(n4s,1),1); % case 1: all zeros input
    case 2
        vCR = ones(size(n4s,1),1); % case 2: all ones input
    case 3
        x = (1:25)';
        vCR = computeCRfromP1(n4s,x); % case 3: a conform CR-function
    case 4
        vCR = zeros(size(n4s,1),1);
        vCR(2) = 1; % case 4: all zeros except for edge number 2, which is set to 1
end



averaging_coefficients = computeJ1vCR(c4n,n4e,n4sDb,vCR);

% to get the expected results, compute the Crouzeix-Raviart function again 
% after computing avCR because in computeAvCR the boundary nodes are
% explicitly set to zero
if ismember(test_case, [2,3])
    vCR = computeCRfromP1(n4s,averaging_coefficients);
end

bubble_coefficients = computeJ2vCR(n4s,vCR,averaging_coefficients);


volume_coefficients = computeJ3vCR(c4n,n4e,vCR,averaging_coefficients, bubble_coefficients);

volume_coefficients(21,:)