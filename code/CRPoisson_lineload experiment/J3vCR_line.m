function val = J3vCR_line(n4s_line,~,Gpts4ref,sides_line,averaging_coefficients,bubble_coefficients)

%% this function expresses J3vCR along a line as a function of barycentric coordinates
% see "calculations_integrateJ3vCR_line.pdf" for derivation of the formulas

lambda2 = Gpts4ref;
lambda1 = 1 - lambda2;

% 𝑎1𝜆1 +𝑎2𝜆2
term1 = averaging_coefficients(n4s_line(:,1)) * lambda1 + ...
    averaging_coefficients(n4s_line(:,2)) * lambda2;

% 6𝑏3𝜆1𝜆2
term2 = 6 * bubble_coefficients(sides_line) * lambda1 * lambda2;

val = term1 + term2;

end