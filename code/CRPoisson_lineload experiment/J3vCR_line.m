function val = J3vCR_line(n4s_line,~,Gpts4ref,sides_line,j1,j2)
%% this function expresses J3vCR along a line as a function of barycentric 
%% coordinates. 
% input:  n4s_line   -   (nr_sides_line,2) matrix containing the pair of
%                        node indeces (as in c4n) that make up each side
%                        on the line
%         Gpts4ref   -   double, Gauss point on the reference line (0,1)
%         sides_line -   (nr_sides_line,1) vector that contains the side
%                        indeces (as in n4s) of each side on the line
%         j1         -   (nr_nodes,k) sparse matrix that contains in each
%                        column the coefficients of J1 applied to a 
%                        CR function. 
%         j2         -   (nr_sides,k) spars matrix that contains in each 
%                        column the coefficients of J2 applied to a CR 
%                        function.
% output: val       -    (nr_sides_line,k) sparse matrix that contains
%                        in each column the values of J3vCR at the Gauss
%                        point on each side of the line.
% see "calculations_integrateJ3vCR_line.pdf" for derivation of the formulas

lambda2 = Gpts4ref;
lambda1 = 1 - lambda2;

% 𝑎1𝜆1 +𝑎2𝜆2
term1 = j1(n4s_line(:,1),:) * lambda1 + ...
    j1(n4s_line(:,2),:) * lambda2;

% 6𝑏3𝜆1𝜆2
term2 = 6 * j2(sides_line,:) * lambda1 * lambda2;

val = term1 + term2;

end