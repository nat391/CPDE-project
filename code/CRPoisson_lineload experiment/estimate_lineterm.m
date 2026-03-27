function lineterm4s = estimate_lineterm(c4n,n4s_line, length4s_line, lineload, degree)
%% this function computes the piecewise lineterm h_E * ∫_E l^2 for every 
%% edge on the line
% input: c4n                nr_nodes x 3 matrix from loadGeometry
%        n4s_line           nr_sides_on_line x 2 matrix giving the two
%                           nodes for each side on the line
%        length4s_line      nr_sides_on_line long vector with the length
%                           of each side
%        lineload           function handle defining the lineload to be
%                           applied along the line
%        degree             degree to use for the quadrature rule used to
%                           integrate the lineload over the line
% output: lineterm4s        integral of the lineload over each side that
%                           covers the line
integrand = @(n4s_line,Gpts4p,Gpts4ref) lineload(Gpts4p).^2;

integral4s_line = integrate(c4n,n4s_line,integrand,degree,length4s_line);

lineterm4s = length4s_line .* integral4s_line;
end