function [integral_value] = integrateJ3vCR_lineload(c4n,n4s_line,sides_line,j1,j2,lineload,degree)
%% This function computes the exact integral of J3vCR along a line
%   input: 
%   c4n - coordinates for node
%   n4s_line - nodes for sides along the line
%   sides_line - [nr_lineSides x 1] sides that lie on the line (indexing according to n4s)
%   j1 - from computeJ1vCR
%   j2 - from compute J2vCR
%   lineload - function by which to multiply J3vCR in the line integral
%   output: integral_value - value of the integral along the line

%% comments
%  1. the volume coefficients are not needed because they are zero along
%  the edges
%  2. we assume that the mesh "matches" the line



% define the integrand J3vCR(lambda1, lambda2, lambda3)
integrand = @(n4s_line,Gpts4p,Gpts4ref) ...
            J3vCR_line(n4s_line,Gpts4p,Gpts4ref,sides_line,j1,j2) ...
            .* lineload(Gpts4p);

integral_value = sum(integrate(c4n,n4s_line,integrand,degree), 'all');