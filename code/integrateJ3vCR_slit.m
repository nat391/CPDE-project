function [integral_value] = integrateJ3vCR_slit(c4n,n4s,averaging_coefficients,bubble_coefficients,s4slit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% assumption: The mesh "matches" the slit


% find nodes and sides that are on the slit
% mask = find(c4n(:,2) == 0);
% nodes_slit = find(c4n(mask,1) >= 0);
% sides_slit = zeros(2^level, 1);
% nr_vertices = size(c4n,1);
% for i = 1:nr_vertices
%     sides_slit(i) = s4n(nodes_slit(i), nodes(i+1))
% 
% end

% compute length4s_slit
nr_slitSides = size(s4slit,1);
length4s_slit = zeros(nr_slitSides,1);
for i = 1:nr_slitSides
    curNode1 = n4s(s4slit,1);
    curNode2 = n4s(s4slit,2);
    length4s_slit(i) = norm(c4n(curNode1,:)' - c4n(curNode2,:)'); %apply norm row-wise on the matrix or write it in a vecotrized form to avoid for loop?
end

integral_value_pw = length4s_slit .* ...
    (sum(averaging_coefficients(n4s(s4slit)),2) / 2 ...
    + bubble_coefficients(s4slit));

integral_value = sum(integral_value_pw);