function [integral_value] = integrateJ3vCR(averaging_coefficients, bubble_coefficients, volume_coefficients, c4n, n4e)
%% function to compute the integral of J3vCR over the whole triangulation
%  Input:  c4n,n4e                  - mesh
%          averaging_coefficients   - coeffs from J1; size: nr_vertices x 1
%          bubble_coefficients      - coeffs from J2; size: nr_sides x 1
%          volume_coefficients      - coeffs from J3; size: nr_elements x 3
%
%  Output: integral_value     - the exact integral of J3vCR over the domain

% compute sides4elements and area4elements matrix
s4e = computeS4e(n4e);
area4e = computeArea4e(c4n, n4e);

% define constants alpha and beta
alpha = sqrt(20 / 27) * (sqrt(7) + 1);
beta = sqrt(20 / 27) * 3 * sqrt(7);

% compute the local sum of the coefficients on each triangle
averaging_coefficients_sum = sum(averaging_coefficients(n4e), 2);
bubble_coefficients_sum = sum(bubble_coefficients(s4e), 2);
volume_coefficients_sum = sum(volume_coefficients, 2);

% compute integral on every triangle
piecewise_integral = area4e .* ... 
    ( averaging_coefficients_sum / 3 + bubble_coefficients_sum / 2 ) ...
    + sqrt(area4e) * ( 9 * alpha - 3 * beta ) .* volume_coefficients_sum / 20;

% sum up the piecewise integrals
integral_value = sum(piecewise_integral,1);

end