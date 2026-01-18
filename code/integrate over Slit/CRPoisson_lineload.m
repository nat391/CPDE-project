function CRPoisson_lineload(max_level, lineload)
% analysing the Poisson-model problem with a line-load with J3vCR

%% initialization

% initialize mesh and useful helpers
[c4n, n4e, n4sDb, n4sNb] = loadGeometry('BigSquare');
n4s = computeN4s(n4e);
s4e = computeS4e(n4e);
[p4n, counts_per_node] = computeP4n(n4e,size(c4n,1));

% find nodes and sides on the slit
nodes_line = find(c4n(:,1)>=0 & c4n(:,2)==0);
sides_line = find( all( ismember(n4s, nodes_line), 2 ) );
n4s_line = n4s(sides_line,:);

% get degrees of freedom -> non-dirichlet sides
n4s = sort(n4s,2);
n4sDb = sort(n4sDb,2);
[isMatch, ~] = ismember(n4s,n4sDb,'rows');
dofs = find(~isMatch);

% initialize error4lvl and nrDof4lvl
error4lvl = zeros(max_level,1);
nrDof4lvl = zeros(max_level,1);

% set Dirichlet and Neumann boundary conditions
u4Db = @(x) 0;
g = @(x) 0;
f = @(x) 1;

for level = 1:max_level

    % useful constants
    nr_sides = size(n4s,1);
    nr_elements = size(n4e,1);

    %% solve

    % compute the RHS
    b = zeros(nr_sides,1);

    for j = dofs'

        %define vCR = psi_j
        vCR = zeros(nr_sides,1);
        vCR(j) = 1;
        
        %compute coefficients for J3vCR    
        averaging_coefficients = computeJ1vCR(s4e,p4n,n4sDb,counts_per_node,vCR);
        bubble_coefficients = computeJ2vCR(n4s,vCR,averaging_coefficients);
        
        % intergrate J3vCR over the slit
        b(j) = integrateJ3vCR_lineload(c4n,n4s_line,sides_line,averaging_coefficients,bubble_coefficients,lineload);
    end
    
    [x, nrDof4lvl(level)] = solveCRPoisson_exactRHS(b,g,u4Db,c4n,n4e,n4sDb,n4sNb);

    %% estimate
    [eta4s, ~] = estimateCREtaSides_noNeumann(f,g,u4Db,x,c4n,n4e,n4sDb);
    error4lvl(level) = sqrt(sum(eta4s));
    % grad4e  = zeros(nr_elements,2);
    % for elem = 1 : nr_elements
    %     grads = 4 * ([1,1,1;c4n(n4e(elem,:),:)'] \ [0,0;eye(2)]);
    %     grad4e(elem,:) = x(s4e(elem,:))'*grads;
    % end
    % eta4e = estimateSigmaAveragingP1(c4n,n4e,grad4e);
    %error4lvl(level) = sqrt(sum(eta4e));
    error4lvl(level) = sqrt(sum(eta4s));

    if(level < max_level)
        %% refine
        [c4n,n4e,n4sDb,n4sNb,n4s_line] = refineUniformRed_slit(c4n,n4e,n4s,n4sDb,n4sNb,n4s_line);
        n4s = computeN4s(n4e); %optimization idea: compute n4s in refinement step
        s4e = computeS4e(n4e);
        [p4n,counts_per_node] = computeP4n(n4e,size(c4n,1));
        sides_line = computeSides_line(n4s,n4s_line);
        

        % get new degrees of freedom -> non-dirichlet sides
        n4s = sort(n4s,2);
        n4sDb = sort(n4sDb,2);
        [isMatch, ~] = ismember(n4s,n4sDb,'rows');
        dofs = find(~isMatch);
    end

end

% plot convergence
figure;
plotConvergence(nrDof4lvl, error4lvl, "F(J3vCR)")
hold on
loglog(nrDof4lvl,nrDof4lvl.^(-0.333))
figure;
plotCR(c4n,n4e,x,'CR-solution')


