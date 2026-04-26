function val = uExact_1(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% get x and y coordinates as row vectors
x1 = x(:,1)';
x2 = x(:,2)';

nr_terms = 50;
n = (1:nr_terms)';

factor1 = (2 * (cos((n*pi)/2) - (-1).^n)) ./ ((n*pi).^2 .* cosh((n*pi)/2));

factor2 = sinh( (n*pi)/2 * (1-abs(x2)) );

factor3 = sin( (n*pi)/2 * (x1+1));

val_matrix = factor1 .* factor2 .* factor3;

% sum over n=1,...,nr_terms and convert back to a column vector
val = sum(val_matrix,1)';

end