function [sides_line] = computeSides_line(n4s,n4s_line)
% This function computes the sides on the line described be the nodes in 
% n4s_line. 
% input:  n4s           nr_sides x 3 matrix from computeN4s
%         n4s_line      nr_sides_on_line x 2 matrix giving the two
%                           nodes for each side on the line
%
% output: sides_line   nr_sides_on_line long vector with the side numbers
%                      for the sides on the line
n4s = sort(n4s,2);
n4s_line = sort(n4s_line,2);
sides_line = find(ismember(n4s,n4s_line,'rows'));
end