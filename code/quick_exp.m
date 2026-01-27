[c4n, n4e, n4sDb, n4sNb] = loadGeometry('BigSquare');

max_level = 9;
ndof4lvl = zeros(max_level);
nr_triangles = zeros(max_level);
nr_edges = zeros(max_level);

for level = 1:max_level
    n4s = computeN4s(n4e);
    n4s = sort(n4s,2);
    n4sDb = sort(n4sDb,2);
    [isMatch, ~] = ismember(n4s,n4sDb,'rows');
    dofs = find(~isMatch);

    ndof4lvl(level) = size(dofs,1);
    nr_triangles(level) = 8*4^(level-1);
    nr_edges(level) = 2^(level-1);

    [c4n,n4e,n4sDb,n4sNb] = refineUniformRed(c4n,n4e,n4sDb,n4sNb);
end



figure;
loglog(ndof4lvl, nr_triangles, 'r-');   % red
hold on
loglog(ndof4lvl, nr_edges, 'b-');       % blue
hold off
grid on
