function main
max_level = 7;

% first experiment
lineload = @(x) 1;
degree = 0;
CRPoisson_lineload_old(max_level,lineload,degree, "l_1");

% second experiment
lineload = @(x) 1-2*x(:,1);
degree = 2;
CRPoisson_lineload_old(max_level,lineload,degree,"l_2" );

%third experiment
lineload = @(x) exp(-x(:,1));
degree = 10;
CRPoisson_lineload_old(max_level,lineload,degree,"l_3");

end