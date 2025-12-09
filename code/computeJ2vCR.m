function [bubble_coefficients] = computeJ2vCR(n4s,vCR,averaging_coefficients)
% this function computes the coefficients J2vCR for the J2 smoothing of a
% Crouzeix-Raviart function vCR
% input: n4s matrix, vCR Crouzeix-Raviart coefficient vector, 
% averaging_coefficients vector of the J1 smoothing of vCR

% compute integral mean of (vCR - avCR) over edge F using midpoint rule
bubble_coefficients = vCR - sum(averaging_coefficients(n4s),2) / 2;

end