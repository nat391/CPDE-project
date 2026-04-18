max_level = 7;
geometry = 'BigSquare';
c4n = loadGeometry(geometry);
nodes_line = find(c4n(:,1)>=0 & c4n(:,2)==0);


ndof4lvl = zeros(max_level,1);
nsides_nodepatch4lvl = zeros(max_level,1);

% initialize mesh and useful helpers
[c4n, n4e, n4sDb, ~] = loadGeometry(geometry);
nr_nodes = size(c4n,1);
n4s = computeN4s(n4e);
s4e = computeS4e(n4e);
[pos4n, counts_per_node] = computePos4n(n4e,nr_nodes);


% find nodes and sides on the line and get their length
sides_line = find( all( ismember(n4s, nodes_line), 2 ) );
n4s_line = n4s(sides_line,:);
length4s_line = computeLength4s(c4n,n4s_line);

for level = 1:max_level

    % get sides of the nodepatches around the line nodes
    sides_nodepatch_line = computeSides_nodepatch(nodes_line,pos4n,s4e);
    nsides_nodepatch4lvl(level) = size(sides_nodepatch_line,1);

    % get degrees of freedom
    n4s = sort(n4s,2);
    n4sDb = sort(n4sDb,2);
    [isMatch, ~] = ismember(n4s,n4sDb,'rows');
    dofs = find(~isMatch);
    ndof4lvl(level) = size(dofs,1);

    %% refine
    [c4n,n4e,n4sDb,n4s_line] = refineUniformRed_with_line(c4n,n4e,n4s,n4sDb,n4s_line);
    nr_nodes = size(c4n,1);
    n4s = computeN4s(n4e); %optimization idea: compute n4s in refinement step
    s4e = computeS4e(n4e);
    [pos4n,counts_per_node] = computePos4n(n4e,nr_nodes);
    nodes_line = unique(n4s_line(:)); % not optimal
    sides_line = computeSides_line(n4s,n4s_line);
    length4s_line = computeLength4s(c4n,n4s_line);
end

% Assume nsides4lvl and ndof4lvl are already defined (1x8 or 8x1)
x = 1:8;

figure;
plot(x, nsides_nodepatch4lvl, '-o', 'LineWidth', 1.5);
hold on;
plot(x, ndof4lvl, '-s', 'LineWidth', 1.5);
hold off;

xlabel('Level');
ylabel('Value');
title('nsides_nodepatch4lvl and ndof4lvl over Levels 1–8');
legend('nsides_nodepatch4lvl', 'ndof4lvl', 'Location', 'best');
grid on;