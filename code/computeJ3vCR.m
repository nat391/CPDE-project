function [j3] = computeJ3vCR(n4e,s4e,area4e,vCR, j1, j2)
% This function applies the J3 operator to a given Crouzeix-Raviart
% function vCR. 
% input: c4n                nr_nodes x 3 matrix from loadGeometry
%        n4e                nr_elems x 3 matrix from loadGeometry
%        vCR                nr_sides long vector with coefficients 
%                           of the CR-function at the edge-midpoints
%        j1                 nr_nodes long vector from computeJ1vCR
%        j2                 nr_sides long vector from computeJ2vCR
%
% Output: j3                nr_elems x 3 matrix with the coefficients of
%                           J3vCR for the three nodes of every triangle

% reorder matrix such that side j is opposite to node j 
% (local enumeration with j=1,2,3)
s4e = s4e(:,[2 3 1]);

% define constants alpha and beta
alpha = sqrt(20 / 27) * (sqrt(7) + 1);
beta = sqrt(20 / 3) * sqrt(7);

% sum up the coefficients for each triangle and copy into three columns
vCR_sum = sum(vCR(s4e),2) * ones(1,3);
averaging_coefficients_sum = sum(j1(n4e),2) * ones(1,3);
bubble_coefficients_sum = sum(j2(s4e),2) * ones(1,3);

% define helper function
computeOther = @(A) [ ...
    sum(A(:,[2,3]),2), ...
    sum(A(:,[3,1]),2), ...
    sum(A(:,[1,2]),2) ...
];

% for each node/edge, sum over the coefficients on the two other
% nodes/edges
vCR_other = computeOther(vCR(s4e));
averaging_coefficients_other = computeOther(j1(n4e));
bubble_coefficients_other = computeOther(j2(s4e));


% terms in the integral
P0_and_P1_term = vCR_sum * (alpha - beta) - alpha * averaging_coefficients_sum;
P2_case1_term = beta * (2*vCR_other + averaging_coefficients_other) - 6 * alpha * bubble_coefficients_sum;
P2_case2_term = beta * (2 * vCR(s4e) + j1(n4e));
P3_case1_term = beta * j2(s4e);
P3_case2_term = beta * bubble_coefficients_other;

% using integration formula for barycentric coordinates
% The derivation of the formula is outlined in
%   "calculations_J3operator.pdf"
j3 = sqrt(area4e) .* ... 
    (P0_and_P1_term / 3 ...
     + P2_case1_term / 12 ...
     + P2_case2_term / 6 ...
     + P3_case1_term / 10 ...
     + P3_case2_term / 5 ...
     );

end