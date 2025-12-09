function [vCR] = computeCRfromP1(n4s,x)
% this function computes the Crouzeix-Raviart coefficients from a given P1
% function


vCR = sum(x(n4s),2) / 2; %calculate the average over each pair of nodes

end