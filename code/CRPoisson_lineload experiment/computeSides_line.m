function [sides_line] = computeSides_line(n4s,n4s_line)
n4s = sort(n4s,2);
n4s_line = sort(n4s_line,2);
sides_line = find(ismember(n4s,n4s_line,'rows'));
end