%% Script, comparing different rhs F(psi_j) and F(J_3 psi_j)
%% in the Poisson-Model-Problem

% import AFEM library
addpath(genpath('C:\Users\natha\Documents\MATLAB\afem-base-master'))

% choose shape and maximum refinement-level of the triangulation of Omega
shape = 'Square';
max_level = 6;

% initialize error4lvl and nrDof4lvl
error4lvl = zeros(max_level,1);
nrDof4lvl = zeros(max_level,1);
error4lvl_J3 = zeros(max_level,1);
nrDof4lvl_J3 = zeros(max_level,1);
[c4n, n4e, n4sDb] = loadGeometry(shape, 1);
n4sNb = [];

for level = 1:max_level

    % refine
    if level >= 2
        [c4n, n4e, n4sDb] = refineUniformRed(c4n,n4e,n4sDb,n4sNb);
    end

    % compute useful matrices and constants
    n4s = computeN4s(n4e);
    nr_vertices = size(c4n,1);
    nr_sides = size(n4s,1);
    nr_elements = size(n4e,1);
    
    %% compute rhs
    % initialize rhs-vector tilde(b)
    b = zeros(nr_sides,1);
    for j = setdiff(1:nr_sides, n4sDb)
        %define vCR = psi_j
        vCR = zeros(nr_sides,1);
        vCR(j) = 1;
        
        %compute coefficients for J3vCR    
        % (maybe compute all at once instead of looping over levels)
        averaging_coefficients = computeJ1vCR(c4n, n4e, n4sDb, vCR);
        bubble_coefficients = computeJ2vCR(n4s, vCR, averaging_coefficients);
        volume_coefficients = computeJ3vCR(c4n, n4e, vCR, averaging_coefficients, bubble_coefficients);
        
        % intergrate J3vCR
        b(j) = integrateJ3vCR(averaging_coefficients, bubble_coefficients, volume_coefficients, c4n, n4e);
    end
    
    % solve Poisson-model-problem for rhs b and tilde(b)
    f = @(x) 1;
    g = @(x) 0;
    u4Db = @(x) 0;
    
    [x, nrDof4lvl(level)] = solveCRPoisson(f,g,u4Db,c4n,n4e,n4sDb,n4sNb);
    [x_J3] = solveCRPoisson_exactRHS(b,g,u4Db,c4n,n4e,n4sDb,n4sNb);
    
    % estimate errors
    [eta4s,~] = estimateCREtaSides_noNeummann(f,g,u4Db,x,c4n,n4e,n4sDb);
    [eta4s_J3,~] = estimateCREtaSides_noNeumann(f,g,u4Db,x_J3,c4n,n4e,n4sDb);

    error4lvl(level) = sqrt(sum(eta4s));
    error4lvl_J3(level) = sqrt(sum(eta4s_J3));
end

% plot convergence
figure;
plotConvergence(nrDof4lvl, error4lvl, "F(vCR)")
figure;
plotConvergence(nrDof4lvl, error4lvl_J3, "F(J3 vCR)")