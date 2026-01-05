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
s4n = computeS4n(n4e);
s4e = computeS4e(n4e);
% Dirichlet boundary sides
DbSides = zeros(1,size(n4sDb,1));
for i = 1:size(n4sDb,1)
        DbSides(i) = s4n(n4sDb(i,1),n4sDb(i,2));
end

% initialize error4lvl and nrDof4lvl
error4lvl = zeros(max_level,1);
nrDof4lvl = zeros(max_level,1);

% set Dirichlet and Neumann boundary conditions
u4Db = @(x) 0;
f = @(x) 0; %not used
g = @(x) 0;

for level = 1:max_level

    % useful constants
    nr_vertices = size(c4n,1);
    nr_sides = size(n4s,1);
    nr_elements = size(n4e,1);

    %% compute rhs
    % initialize rhs-vector b
    b = zeros(nr_sides,1);
    for j = setdiff(1:nr_sides, DbSides)
        %define vCR = psi_j
        vCR = zeros(nr_sides,1);
        vCR(j) = 1;
        
        %compute coefficients for J3vCR    
        % (maybe compute all at once instead of looping over levels)
        averaging_coefficients = computeJ1vCR(c4n,n4e, n4sDb, vCR);
        bubble_coefficients = computeJ2vCR(n4s, vCR, averaging_coefficients);
        %bubble_coefficients = zeros(size(bubble_coefficients,1),1);
        
        % intergrate J3vCR
        b(j) = integrateJ3vCR_slit(c4n,n4s,averaging_coefficients,bubble_coefficients,s4slit);
    end
    %solve
    [x, nrDof4lvl(level)] = solveCRPoisson_exactRHS(b,g,u4Db,c4n,n4e,n4sDb,n4sNb);

    % estimate
    %[eta4s, ~] = estimateCREtaSides_noNeumann(f,g,u4Db,x,c4n,n4e,n4sDb);
    %error4lvl(level) = sqrt(sum(eta4s));
    grad4e  = zeros(nr_elements,2);
    for elem = 1 : nr_elements
        grads = 4 * ([1,1,1;c4n(n4e(elem,:),:)'] \ [0,0;eye(2)]);
        grad4e(elem,:) = x(s4e(elem,:))'*grads;
    end
    eta4e = estimateSigmaAveragingP1(c4n,n4e,grad4e);
    error4lvl(level) = sqrt(sum(eta4e));

    if(level < max_level)
        % refine
        [c4n,n4e,n4sDb,n4sNb,n4sSlit] = refineUniformRed_slit(c4n,n4e,n4s,n4sDb,n4sNb,n4sSlit);
        n4s = computeN4s(n4e); %optimization idea: compute n4s in refinement step
        n4s = sort(n4s,2);
        n4sSlit = sort(n4sSlit,2);
        s4slit = find(ismember(n4s,n4sSlit,'rows'));
        s4n = computeS4n(n4e);
        s4e = computeS4e(n4e);
        % Dirichlet boundary sides
        DbSides = zeros(1,size(n4sDb,1));
        for i = 1:size(n4sDb,1)
                DbSides(i) = s4n(n4sDb(i,1),n4sDb(i,2));
        end
    end

end

% plot convergence
close all
figure;
plotConvergence(nrDof4lvl, error4lvl, "F(J3vCR)")
hold on
loglog(nrDof4lvl,nrDof4lvl.^(-1))
figure;
plotCR(c4n,n4e,x,'CR-solution')


