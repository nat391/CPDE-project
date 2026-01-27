function lineterm4s = estimate_lineterm(c4n,n4s_line, length4s_line, lineload, degree)
%% this function computes the piecewise lineterm h_E * ∫_E l^2 for every 
%% edge on the line
%   Detailed explanation goes here
integrand = @(n4s_line,Gpts4p,Gpts4ref) lineload(Gpts4p).^2;

integral4s_line = integrate(c4n,n4s_line,integrand,degree,length4s_line);

lineterm4s = length4s_line .* integral4s_line;
end