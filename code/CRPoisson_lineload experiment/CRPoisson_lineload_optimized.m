function CRPoisson_lineload_optimized(geometry, nodes_line, max_level, lineload, degree, OPTname)
%% This function solves the Poisson-model problem for a lineload on the rhs
%% using the smoother J. It also plots the convergence rates. 
% input: geometry               string - geometry to use, must be a name of a geom
%                               from AFEMs geometries folder
%        nodes_line             vector of length nr_nodes_on_line with
%                               nodes on the line where lineload is applied
%        max_level              integer - maximum refinement level to comp.
%        lineload               function handle for the lineload that is
%                               applied along the edges connected by
%                               nodes_line
%        degree                 integer - degree to use for the quadrature
%                               rules when integrating lineload
%
% output: none

%% initialization

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

% get sides of the nodepatches around the line nodes
sides_nodepatch_line = computeSides_nodepatch(nodes_line,pos4n,s4e);

% exclude boundary sides
n4s_nodepatch_line = sort(n4s(sides_nodepatch_line,:),2);
n4sDb = sort(n4sDb,2);
[isMatch, ~] = ismember(n4s_nodepatch_line,n4sDb,'rows');
sides_nodepatch_line = sides_nodepatch_line(~isMatch);

% initialize error4lvl and ndof4lvl
error4lvl = zeros(max_level,1);
error4lvl_jumpterm = zeros(max_level,1);
ndof4lvl = zeros(max_level,1);

% set Dirichlet boundary conditions
u4Db = @(x) 0;

for level = 1:max_level

    % useful constants
    nr_sides = size(n4s,1);
    nr_sides_nodepatch_line = length(sides_nodepatch_line);

    %% solve

    % compute the RHS where non-zero
    b_local = zeros(nr_sides_nodepatch_line,1);

    parfor k = 1:nr_sides_nodepatch_line

        side = sides_nodepatch_line(k);

        %compute coefficients for J3vCR    
        j1 = computeJ1psi(s4e,pos4n,n4sDb,counts_per_node,side);
        j2 = computeJ2psi(n4s,side,j1);
        
        % intergrate J3vCR over the slit
        b_local(k) = integrateJ3vCR_lineload(c4n,n4s_line,sides_line,j1,j2,lineload,degree+2);
    end
    
    % scatter back to full rhs vector b
    b = zeros(nr_sides,1);
    b(sides_nodepatch_line) = b_local;
  
    [uCR, ndof4lvl(level)] = solveCRPoisson_exactRHS(b,u4Db,c4n,n4e,n4sDb);

    %% estimate
    [eta4s_jumpterm, ~] = estimate_CRjumpterm(u4Db,uCR,c4n,n4e,n4sDb);
    eta4s_lineterm = estimate_lineterm(c4n,n4s_line,length4s_line,lineload,degree);
    error4lvl_jumpterm(level) = sqrt(sum(eta4s_jumpterm));
    error4lvl(level) = sqrt(sum(eta4s_jumpterm) + sum(eta4s_lineterm));

    if(level < max_level)
        %% refine
        [c4n,n4e,n4sDb,n4s_line] = refineUniformRed_with_line(c4n,n4e,n4s,n4sDb,n4s_line);
        nr_nodes = size(c4n,1);
        n4s = computeN4s(n4e); %optimization idea: compute n4s in refinement step
        s4e = computeS4e(n4e);
        [pos4n,counts_per_node] = computePos4n(n4e,nr_nodes);
        nodes_line = unique(n4s_line(:)); % not optimal
        sides_line = computeSides_line(n4s,n4s_line);
        length4s_line = computeLength4s(c4n,n4s_line);
        
        % get sides of the nodepatches around the line nodes
        sides_nodepatch_line = computeSides_nodepatch(nodes_line,pos4n,s4e);
        
        % exclude boundary sides
        n4s_nodepatch_line = sort(n4s(sides_nodepatch_line,:),2);
        n4sDb = sort(n4sDb,2);
        [isMatch, ~] = ismember(n4s_nodepatch_line,n4sDb,'rows');
        sides_nodepatch_line = sides_nodepatch_line(~isMatch);
        
    end

end

% ---- plot convergence ----
figConv = figure;

plotConvergence(ndof4lvl, error4lvl, "\eta_l");
hold on

plotConvergence(ndof4lvl, error4lvl_jumpterm, "\eta_{l,jump}");

scaling1 = error4lvl(max_level) / (ndof4lvl(max_level)^(-0.25)) * 1.5;
scaling2 = error4lvl_jumpterm(max_level) / (ndof4lvl(max_level)^(-0.5)) * 1.4;

loglog(ndof4lvl, scaling1 * ndof4lvl.^(-0.25), 'k--', 'DisplayName', 's=-0.25');
loglog(ndof4lvl, scaling2 * ndof4lvl.^(-0.5), 'k--', 'DisplayName', 's=-0.5');

xlabel('ndof')
legend('show');


% ---- plot solution ----
figSol = figure;
plotCR(c4n,n4e,uCR,OPTname);

% ---- save figures ----
% outDir = "images";
% if ~exist(outDir, "dir")
%     mkdir(outDir);
% end





