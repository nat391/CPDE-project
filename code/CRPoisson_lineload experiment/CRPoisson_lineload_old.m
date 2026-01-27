function CRPoisson_lineload_old(max_level, lineload, degree, OPTname)
%% This function solves the Poisson-model problem for a lineload on the rhs
%% using the smoother J3. It also plots the convergence rates. 
% It uses the BigSquare domain and the line [0,1] x {0}. 

if nargin < 4 || isempty(OPTname)
        baseName = 'Figure';   % default name
else
        baseName = OPTname;
end


%% initialization

% initialize mesh and useful helpers
[c4n, n4e, n4sDb, n4sNb] = loadGeometry('BigSquare');
nr_nodes = size(c4n,1);
n4s = computeN4s(n4e);
s4e = computeS4e(n4e);
[pos4n, counts_per_node] = computePos4n(n4e,nr_nodes);


% find nodes and sides on the line and get their length
nodes_line = find(c4n(:,1)>=0 & c4n(:,2)==0);
sides_line = find( all( ismember(n4s, nodes_line), 2 ) );
n4s_line = n4s(sides_line,:);
length4s_line = computeLength4s(c4n,n4s_line);

% get degrees of freedom -> big patch sides without boundary sides
n4s = sort(n4s,2);
n4sDb = sort(n4sDb,2);
[isMatch, ~] = ismember(n4s,n4sDb,'rows');
dofs = find(~isMatch);

% initialize error4lvl and ndof4lvl
error4lvl = zeros(max_level,1);
ndof4lvl = zeros(max_level,1);

% set Dirichlet and Neumann boundary conditions
u4Db = @(x) 0;
g = @(x) 0;

for level = 1:max_level

    % useful constants
    nr_sides = size(n4s,1);

    %% solve

    % compute the RHS
    b = zeros(nr_sides,1);

    for j = dofs'

        %define vCR = psi_j
        vCR = zeros(nr_sides,1);
        vCR(j) = 1;
        
        %compute coefficients for J3vCR    
        averaging_coefficients = computeJ1vCR(s4e,pos4n,n4sDb,counts_per_node,vCR);
        bubble_coefficients = computeJ2vCR(n4s,vCR,averaging_coefficients);
        
        % intergrate J3vCR over the slit
        b(j) = integrateJ3vCR_lineload(c4n,n4s_line,sides_line,averaging_coefficients,bubble_coefficients,lineload,degree+2);
    end
    
    [x, ndof4lvl(level)] = solveCRPoisson_exactRHS(b,g,u4Db,c4n,n4e,n4sDb,n4sNb);

    %% estimate
    [eta4s_jumpterm, ~] = estimateCREtaSides_noNeumann(g,u4Db,x,c4n,n4e,n4sDb);
    eta4s_lineterm = estimate_lineterm(c4n,n4s_line,length4s_line,lineload,degree);
    error4lvl(level) = sqrt(sum(eta4s_jumpterm) + sum(eta4s_lineterm));

    if(level < max_level)
        %% refine
        [c4n,n4e,n4sDb,n4sNb,n4s_line] = refineUniformRed_with_line(c4n,n4e,n4s,n4sDb,n4sNb,n4s_line);
        nr_nodes = size(c4n,1);
        n4s = computeN4s(n4e); %optimization idea: compute n4s in refinement step
        s4e = computeS4e(n4e);
        [pos4n,counts_per_node] = computePos4n(n4e,nr_nodes);
        sides_line = computeSides_line(n4s,n4s_line);
        length4s_line = computeLength4s(c4n,n4s_line);
        

        % get new degrees of freedom
        n4s = sort(n4s,2);
        n4sDb = sort(n4sDb,2);
        [isMatch, ~] = ismember(n4s,n4sDb,'rows');
        dofs = find(~isMatch);
    end

end

% ---- plot convergence ----
figConv = figure;
plotConvergence(ndof4lvl, error4lvl, "\eta_l");

% ---- plot solution ----
figSol = figure;
plotCR(c4n,n4e,x,'CR-solution');

% ---- save figures ----
outDir = "images";
if ~exist(outDir, "dir")
    mkdir(outDir);
end

%exportgraphics(figConv, fullfile(outDir, baseName + "_convergence.png"), "Resolution", 300);
%exportgraphics(figSol,  fullfile(outDir, baseName + "_solution.png"),     "Resolution", 300);

% optional: also save .fig for later editing in MATLAB
%savefig(figConv, fullfile(outDir, baseName + "_convergence.fig"));
%savefig(figSol,  fullfile(outDir, baseName + "_solution.fig"));





