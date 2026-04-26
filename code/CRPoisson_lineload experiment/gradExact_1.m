function val = gradExact_1(x)
%% computes the gradient of the exact solution u to the following
%% poisson problem
% Seek for a solution u such that
%   -div(grad(u)) = f   in Omega,
%              u  = 0   on boundary of Omega,
% where Omega is the big square (-1,1)^2 and f(x,y) = theta(x)*delta(y).
% Theta refers to the Heaviside function and delta to the dirac delta
% distribution. 
% Input:  x      -       (nr_points,2) matrix which contains the coordinates
%                        of the points in Omega at which the gradient
%                        should be computed
% Output: val    -       (nr_points,2) matrix, giving the gradient at each
%                        point in x. 

% get x and y coordinates as row vectors
x1 = x(:,1)';
x2 = x(:,2)';

% set the number of terms
nr_terms = 50;
n = (1:nr_terms)';

factor1 = ((cos((n*pi)/2) - (-1).^n)) ./ ((n*pi) .* cosh((n*pi)/2));

%% partial derivative with respect to x1
factor2 = sinh( (n*pi)/2 * (1-abs(x2)) );

factor3 = cos( (n*pi)/2 * (x1+1));

val_matrix = factor1 .* factor2 .* factor3;

% sum over n=1,...,nr_terms and convert back to a column vector
val_partial_x1 = sum(val_matrix,1)';

%% partial derivative with respect to x2
factor2 = -sign(x2) .* cosh( (n*pi)/2 * (1-abs(x2)) );

factor3 = sin( (n*pi)/2 * (x1+1));

val_matrix = factor1 .* factor2 .* factor3;

% sum over n=1,...,nr_terms and convert back to a column vector
val_partial_x2 = sum(val_matrix,1)';

val = [val_partial_x1,val_partial_x2];

end