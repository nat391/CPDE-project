% analysing the integration over a slit with J3vCR

% import AFEM library
addpath(genpath('C:\Users\natha\Documents\MATLAB\afem-base-master'))

% set maximum refinement level
max_level = 5;

% initialize mesh
[c4n, n4e, n4sDb, n4sNb] = loadGeometry('BigSquare');
s4slit = 5;
n4s = computeN4s(n4e);
n4sSlit = n4s(s4slit,:);

% initialize error4lvl and nrDof4lvl
error4lvl = zeros(max_level,1);
nrDof4lvl = zeros(max_level,1);

% set Dirichlet and Neumann boundary conditions
u4Db = @(x) 0;
f = @(x) 0;
g = @(x) 0;

for level = 1:max_level

    % useful constants
    nr_vertices = size(c4n,1);
    nr_sides = size(n4s,1);
    nr_elements = size(n4e,1);

    %% compute rhs
    % initialize rhs-vector b
    b = zeros(nr_sides,1);
    for j = setdiff(1:nr_sides, n4sDb)
        %define vCR = psi_j
        vCR = zeros(nr_sides,1);
        vCR(j) = 1;
        
        %compute coefficients for J3vCR    
        % (maybe compute all at once instead of looping over levels)
        averaging_coefficients = computeJ1vCR(c4n,n4e, n4sDb, vCR);
        bubble_coefficients = computeJ2vCR(n4s, vCR, averaging_coefficients);
        
        % intergrate J3vCR
        b(j) = integrateJ3vCR_slit(c4n,n4s,averaging_coefficients,bubble_coefficients,s4slit);
    end
    %solve
    [x, nrDof4lvl(level)] = solveCRPoisson_exactRHS(b,g,u4Db,c4n,n4e,n4sDb,n4sNb);

    % estimate
    [eta4s, ~] = estimateCREtaSides_noNeumann(f,g,u4Db,x,c4n,n4e,n4sDb);
    error4lvl(level) = sqrt(sum(eta4s));

    % refine
    [c4n,n4e,n4sDb,n4sNb,n4sSlit] = refineUniformRed_slit(c4n,n4e,n4s,n4sDb,n4sNb,n4sSlit);
    n4s = computeN4s(n4e); %optimization idea: compute n4s in refinement step
    n4s = sort(n4s,2);
    n4sSlit = sort(n4sSlit,2);
    s4slit = find(ismember(n4s,n4sSlit,'rows'));
end

% plot convergence
figure;
plotConvergence(nrDof4lvl, error4lvl, "F(vCR)")


