function main

%% This function runs the experiment described in the project report
% The domain is the big square (-1,1)^2 with the line (0,1)x{0}. The three
% lineloads are 1, 1-2x and exp(-x). 
max_level = 10;
geometry = 'BigSquare';
c4n = loadGeometry(geometry);
nodes_line = find(c4n(:,1)>=0 & c4n(:,2)==0);

% first experiment
lineload = @(x) 1;
degree = 0;
CRPoisson_lineload_optimized(geometry,nodes_line,max_level,lineload,degree, "l_1");

% second experiment
lineload = @(x) 1-2*x(:,1);
degree = 2;
CRPoisson_lineload_optimized(geometry,nodes_line,max_level,lineload,degree,"l_2" );

%third experiment
lineload = @(x) exp(-x(:,1));
degree = 10;
CRPoisson_lineload_optimized(geometry,nodes_line,max_level,lineload,degree,"l_3");

end