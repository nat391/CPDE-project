function CRPoisson_lineload_adaptive(geometry, nodes_line, min_ndof, lineload, degree)
%% This function solves the Poisson-model problem for a lineload on the rhs
%% using the smoother J. It also plots the convergence rates. 
% input: geometry               string - geometry to use, must be a name of 
%                               a geometry from AFEMs geometries folder.
%        nodes_line             vector of length nr_nodes_on_line with
%                               nodes on the line where lineload is applied
%        min_ndof               integer - minimum degrees of freedom to
%                               reach before terminating
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
dirichlet_nodes = unique(n4sDb(:));
n4s = computeN4s(n4e);
s4e = computeS4e(n4e);
[pos4n, counts_per_node] = computePos4n(n4e,nr_nodes);


% find nodes and sides on the line and get their length
sides_line = find(all(ismember(n4s,nodes_line),2));
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
eta4ndof = sparse(1,1);
eta4ndof_lineterm = sparse(1,1);
eta4ndof_jumpterm = sparse(1,1);

% set Dirichlet boundary conditions
u4Db = @(x) 0;

while (true)

    % useful constants
    nr_sides = size(n4s,1);
    nr_sides_nodepatch_line = length(sides_nodepatch_line);

    %% solve
    
    b = zeros(nr_sides,1);

    %compute coefficients for J3vCR in a vectorized manner
    j1_all = computeJ1psi(s4e,pos4n,dirichlet_nodes,counts_per_node,...
                                                     sides_nodepatch_line);
    j2_all = computeJ2psi(n4s,sides_nodepatch_line,j1_all);
    b(sides_nodepatch_line) = integrateJ3psi_lineload(...
                  c4n,n4s_line,sides_line,j1_all,j2_all,lineload,degree+2);
  
    [uCR,ndof] = solveCRPoisson_exactRHS(b,u4Db,c4n,n4e,n4sDb);

    %% estimate
    eta4s_lineterm = estimate_lineterm(c4n,n4s_line,length4s_line,...
                                                          lineload,degree);
    eta4s_jumpterm = estimate_CRjumpterm(u4Db,uCR,c4n,n4e,n4s,s4e,n4sDb);
    
    eta4ndof_lineterm(ndof) = sqrt(sum(eta4s_lineterm));
    eta4ndof_jumpterm(ndof) = sqrt(sum(eta4s_jumpterm));
    eta4ndof(ndof) = sqrt(sum(eta4s_jumpterm) + sum(eta4s_lineterm));
    
    % break if minimum degrees of freedom were reached. Else, continue.
    if ndof >= min_ndof, break, end

    %% mark
    eta4s = zeros(nr_sides,1);
    eta4s(sides_line) = eta4s_lineterm;
    eta4s = eta4s + eta4s_jumpterm;
    n4sMarked = markBulk(n4s,eta4s);

    %% refine
    [c4n,n4e,nodes_line,n4sDb] = refineRGB_with_line(c4n,n4e,nodes_line,sides_line,n4sDb,n4sMarked);
    
    % compute helpers again
    nr_nodes = size(c4n,1);
    n4s = computeN4s(n4e);
    s4e = computeS4e(n4e);
    [pos4n,counts_per_node] = computePos4n(n4e,nr_nodes);

    % find nodes and sides on the line and get their length
    sides_line = find(all(ismember(n4s,nodes_line),2));
    n4s_line = n4s(sides_line,:);
    length4s_line = computeLength4s(c4n,n4s_line);
    
    % get sides of the nodepatches around the line nodes
    sides_nodepatch_line = computeSides_nodepatch(nodes_line,pos4n,s4e);
    
    % exclude boundary sides
    n4s_nodepatch_line = sort(n4s(sides_nodepatch_line,:),2);
    n4sDb = sort(n4sDb,2);
    [isMatch, ~] = ismember(n4s_nodepatch_line,n4sDb,'rows');
    sides_nodepatch_line = sides_nodepatch_line(~isMatch);

end

% ---- plot triangulation ----
figTriangulation = figure;
plotTriangulation(c4n,n4e);

% ---- plot solution ----
figSolution = figure;
plotCR(c4n,n4e,uCR);

% ---- plot convergence ----
ndof4lvl = find(eta4ndof);
eta4lvl_lineterm = eta4ndof_lineterm(ndof4lvl);
eta4lvl_jumpterm = eta4ndof_jumpterm(ndof4lvl);
eta4lvl = eta4ndof(ndof4lvl);
figConvergence = figure;

plotConvergence(ndof4lvl, eta4lvl_lineterm, "\eta_{l,line}");
hold on

plotConvergence(ndof4lvl, eta4lvl_jumpterm, "\eta_{l,jump}");
hold on

plotConvergence(ndof4lvl, eta4lvl, "\eta_l");
hold on

%scaling1 = error4lvl(max_level) / (ndof4lvl(max_level)^(-0.25)) * 1.5;
scaling2 = eta4ndof_jumpterm(end) / (ndof4lvl(end)^(-0.5)) * 1.4;

%loglog(ndof4lvl, scaling1 * ndof4lvl.^(-0.25), 'k--', 'DisplayName', 's=-0.25');
loglog(ndof4lvl, scaling2 * ndof4lvl.^(-0.5), 'k--', 'DisplayName', 's=-0.5');

xlabel('ndof')
legend('show');


% ---- save figures ----
% Ordner "plots" relativ zum aktuellen Skript anlegen (falls nicht vorhanden)
plotsDir = fullfile(pwd, 'plots');
if ~exist(plotsDir, 'dir')
    mkdir(plotsDir);
end

% Zeitstempel erzeugen (Dateiname-kompatibel)
timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM-SS');

% ---- Figure: Triangulation ----
filenameTri = fullfile(plotsDir, ['Triangulation_' timestamp]);

savefig(figTriangulation, [filenameTri '.fig']);       % MATLAB Figure
exportgraphics(figTriangulation, [filenameTri '.png'], 'Resolution', 300);

% ---- Figure: Solution ----
filenameSol = fullfile(plotsDir, ['Solution_' timestamp]);

savefig(figSolution, [filenameSol '.fig']);       % MATLAB Figure
exportgraphics(figSolution, [filenameSol '.png'], 'Resolution', 300);

% ---- Figure: Convergence ----
filenameConv = fullfile(plotsDir, ['Convergence_' timestamp]);

savefig(figConvergence, [filenameConv '.fig']);   % MATLAB Figure
exportgraphics(figConvergence, [filenameConv '.png']);

disp('Plots erfolgreich gespeichert.');