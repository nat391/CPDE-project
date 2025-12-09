function [volume_coefficients] = computeJ3vCR(c4n,n4e,vCR, averaging_coefficients, bubble_coefficients)
% This function applies the J3 operator onto a given Crouzeix-Raviart
% function vCR. 
% Input: c4n, n4e, vCR, averaging coefficients and bubble coefficients
% Output: volume coefficients, a nr_elements x 3 matrix
%   The derivation of the formula is outlined in
%   "calculations_J3operator.pdf"

% compute the |cal(T)| x 3 "sides for elements" matrix
s4e = computeS4e(n4e);
% reorder matrix such that side j is opposite to node j 
% (local enumeration with j=1,2,3)
s4e = s4e(:,[2 3 1]);
Area4e = computeArea4e(c4n,n4e);

% define constants alpha and beta
alpha = sqrt(20 / 27) * (sqrt(7) + 1);
beta = sqrt(20 / 27) * 3 * sqrt(7);

% sum up the coefficients for each triangle and copy into three columns
vCR_sum = sum(vCR(s4e),2) * ones(1,3);
averaging_coefficients_sum = sum(averaging_coefficients(n4e),2) * ones(1,3);
bubble_coefficients_sum = sum(bubble_coefficients(s4e),2) * ones(1,3);

% define helper function
computeOther = @(A) [ ...
    sum(A(:,[2,3]),2), ...
    sum(A(:,[3,1]),2), ...
    sum(A(:,[1,2]),2) ...
];

% for each node/edge, sum over the coefficients on the two other
% nodes/edges
vCR_other = computeOther(vCR(s4e));
averaging_coefficients_other = computeOther(averaging_coefficients(n4e));
bubble_coefficients_other = computeOther(bubble_coefficients(s4e));


% terms in the integral
P0_and_P1_term = vCR_sum * (alpha - beta) - alpha * averaging_coefficients_sum;
P2_case1_term = beta * (2*vCR_other + averaging_coefficients_other) - 6 * alpha * bubble_coefficients_sum;
P2_case2_term = beta * (2 * vCR(s4e) + averaging_coefficients(n4e));
P3_case1_term = beta * bubble_coefficients(s4e);
P3_case2_term = beta * bubble_coefficients_other;

% using integration formula for barycentric coordinates
volume_coefficients = sqrt(Area4e) .* ... 
    (P0_and_P1_term / 3 ...
     + P2_case1_term / 12 ...
     + P2_case2_term / 6 ...
     + P3_case1_term / 10 ...
     + P3_case2_term / 5 ...
     );

end