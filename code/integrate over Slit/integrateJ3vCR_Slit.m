function [integral_value] = integrateJ3vCR_slit(c4n,n4s,averaging_coefficients,bubble_coefficients,sides_slit)
%% This function computes the exact integral of J3 vCR along a slit
%   input: 
%   c4n - coordinates for nodes
%   n4s - nodes for sides
%   averaging_coefficients - from computeJ1vCR
%   bubble_coefficients - from compute J2vCR
%   sides_slit - [nr_slitSides x 1] sides that lie on the slit (indexing according to n4s)
%% comments
%  1. the volume coefficients are not needed because they are zero along
%  the edges
%  2. we assume that the mesh "matches" the slit

%% compute length4s_slit - the length for each side on the slit

% initialize lengths4s_slit
nr_slitSides = size(sides_slit,1);
length4s_slit = zeros(nr_slitSides,1);

% compute distance of two nodes for all sides on the slit
% optimization idea -> take as an input
for i = 1:nr_slitSides
    length4s_slit(i) = norm(c4n(n4s(sides_slit(i),1),:) - c4n(n4s(sides_slit(i),2),:)); %apply norm row-wise on the matrix or write it in a vecotrized form to avoid for loop?
end

% compute the piecewise integral (along each side)
% derivation of the formula in: "calculations_integrateJ3vCR_slit.pdf"
integral_value_pw = length4s_slit .* ...
    (sum(averaging_coefficients(n4s(sides_slit)),2) / 2 ...
    + bubble_coefficients(sides_slit));

% sum up the piecewise integrals
integral_value = sum(integral_value_pw);