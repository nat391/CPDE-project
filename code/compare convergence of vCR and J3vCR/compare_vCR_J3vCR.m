function compare_vCR_J3vCR(max_level, shape)
%% This script compares different rhs F(psi_j) and F(J_3 psi_j)
%% in the Poisson-Model-Problem

%% initialize

% initialize mesh
[c4n, n4e, n4sDb] = loadGeometry(shape, 1);
n4sNb = [];

% initialize boundary conditions
f = @(x) 1;
g = @(x) 0;
u4Db = @(x) 0;

% initialize error4lvl and nrDof4lvl
error4lvl = zeros(max_level,1);
nrDof4lvl = zeros(max_level,1);
error4lvl_J3 = zeros(max_level,1);

% compute useful helpers
n4s = computeN4s(n4e);
s4e = computeS4e(n4e);
[p4n, counts_per_node] = computeP4n(n4e,size(c4n,1));
area4e = computeArea4e(c4n,n4e);

% get degrees of freedom -> non-dirichlet sides
n4s = sort(n4s,2);
n4sDb = sort(n4sDb,2);
[isMatch, ~] = ismember(n4s,n4sDb,'rows');
dofs = find(~isMatch);



for level = 1:max_level

    %% refine
    if level >= 2
        % refine mesh
        [c4n, n4e, n4sDb] = refineUniformRed(c4n,n4e,n4sDb,n4sNb);

        % compute useful helpers
        n4s = computeN4s(n4e);
        s4e = computeS4e(n4e);
        [p4n, counts_per_node] = computeP4n(n4e,size(c4n,1));
        area4e = computeArea4e(c4n,n4e);

        % get new degrees of freedom -> non-dirichlet sides
        n4s = sort(n4s,2);
        n4sDb = sort(n4sDb,2);
        [isMatch, ~] = ismember(n4s,n4sDb,'rows');
        dofs = find(~isMatch);
    end

    %% solve

    % initialize rhs-vector b
    nr_sides = size(n4s,1);
    b = zeros(nr_sides,1);

    % loop over all inner sides
    for j = dofs'

        %define vCR = psi_j
        vCR = zeros(nr_sides,1);
        vCR(j) = 1;
        
        %compute coefficients for J3vCR    
        averaging_coefficients = computeJ1vCR(s4e,p4n,n4sDb,counts_per_node,vCR);
        bubble_coefficients = computeJ2vCR(n4s,vCR,averaging_coefficients);
        volume_coefficients = computeJ3vCR(n4e,s4e,area4e,vCR,averaging_coefficients,bubble_coefficients);
        
        % intergrate J3vCR
        b(j) = integrateJ3vCR(s4e,area4e,averaging_coefficients, bubble_coefficients, volume_coefficients, c4n, n4e);
    end

    % solve Poisson-model-problem
    [x, nrDof4lvl(level)] = solveCRPoisson(f,g,u4Db,c4n,n4e,n4sDb,n4sNb);
    [x_J3] = solveCRPoisson_exactRHS(b,g,u4Db,c4n,n4e,n4sDb,n4sNb);
    figure('Position', [100 100 1000 500])

    % plot solutions
    tiledlayout(1,2)
    nexttile
    plotCR(c4n,n4e,x,"normal solution")
    nexttile
    plotCR(c4n,n4e,x_J3,"J3 solution")
    drawnow
    
    %% estimate
    [eta4s,~] = estimateCREtaSides_noNeumann(f,g,u4Db,x,c4n,n4e,n4sDb);
    [eta4s_J3,~] = estimateCREtaSides_noNeumann(f,g,u4Db,x_J3,c4n,n4e,n4sDb);

    error4lvl(level) = sqrt(sum(eta4s));
    error4lvl_J3(level) = sqrt(sum(eta4s_J3));
end

% plot convergence
figure;
plotConvergence(nrDof4lvl, error4lvl, "F(vCR)")
hold on;
plotConvergence(nrDof4lvl, error4lvl_J3, "F(J3 vCR)")
hold on;
loglog(nrDof4lvl, 10*nrDof4lvl.^(-1/2), 'DisplayName', "slope = -1/2")