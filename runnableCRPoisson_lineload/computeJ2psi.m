function [j2] = computeJ2psi(n4s,side,j1)
% this function computes the coefficients j2 of the J2 operator applied 
% to a CR-function vCR. 
% input:  n4s             nr_sides x 2 matrix from computeN4s
%         side            side index where CR basis function is 1 
%         j1              nr_nodes long vector from computeJ1vCR
%
% output: j2             nr_sides long vector with the coefficients 
%                        of degrees of freedom that J2vCR has along the edges

% compute integral mean of (vCR - J1vCR) over edge F using midpoint rule
j2 = -sum(j1(n4s),2) / 2;
j2(side) = j2(side) + 1;

end
