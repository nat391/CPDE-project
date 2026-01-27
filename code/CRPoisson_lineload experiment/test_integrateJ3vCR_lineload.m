function tests = test_integrateJ3vCR_lineload
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    refinement_lvl = 5;
    [c4n,n4e,n4sDb] = loadGeometry('BigSquare',refinement_lvl);
    
    % find nodes that are on the slit
    n4s_line = find(c4n(:,1)>=0 & c4n(:,2)==0);
    n4s = computeN4s(n4e);

    % get the indeces of the sides that are on the slit
    sides_line = find( all( ismember(n4s, n4s_line), 2 ) );
    n4s_line = n4s(sides_line,:);

    % other useful matrices
    s4e = computeS4e(n4e);
    [pos4n, counts_per_node] = computePos4n(n4e,size(c4n,1));

    % store everything
    testCase.TestData.c4n        = c4n;
    testCase.TestData.n4e        = n4e;
    testCase.TestData.n4s       = n4s;
    testCase.TestData.n4sDb     = n4sDb;
    testCase.TestData.sides_line = sides_line;
    testCase.TestData.n4s_line = n4s_line;
    testCase.TestData.s4e = s4e;
    testCase.TestData.pos4n = pos4n;
    testCase.TestData.counts_per_node = counts_per_node;
end


function testZeroInput(testCase)

    % retrieve shared data
    c4n        = testCase.TestData.c4n;
    n4s        = testCase.TestData.n4s;
    n4sDb      = testCase.TestData.n4sDb;
    sides_line = testCase.TestData.sides_line;
    n4s_line = testCase.TestData.n4s_line;
    s4e = testCase.TestData.s4e;
    pos4n = testCase.TestData.pos4n;
    counts_per_node = testCase.TestData.counts_per_node;
    

    expSolution = 0;

    lineload = @(x) 1;
    
    nr_sides = size(n4s,1);
    vCR = zeros(nr_sides,1);

    averaging_coefficients = computeJ1vCR(s4e,pos4n,n4sDb,counts_per_node,vCR);
    bubble_coefficients   = computeJ2vCR(n4s,vCR,averaging_coefficients);

    actSolution = integrateJ3vCR_lineload(...
        c4n,n4s_line,sides_line,averaging_coefficients,bubble_coefficients,lineload);

    verifyEqual(testCase, expSolution, actSolution)
end


function testOneInput(testCase)

    expSolution = 1;

    % define lineload
    lineload = @(x) 1;

    % retrieve shared data
    c4n        = testCase.TestData.c4n;
    n4s        = testCase.TestData.n4s;
    n4sDb      = testCase.TestData.n4sDb;
    sides_line = testCase.TestData.sides_line;
    n4s_line = testCase.TestData.n4s_line;
    s4e = testCase.TestData.s4e;
    pos4n = testCase.TestData.pos4n;
    counts_per_node = testCase.TestData.counts_per_node;

    nr_sides = size(n4s,1);
    
    % find dofs -> sides that are not on the boundary
    n4s = sort(n4s,2);
    n4sDb = sort(n4sDb,2);
    [isMatch, ~] = ismember(n4s,n4sDb,'rows');
    
    % define CR-input
    vCR = zeros(nr_sides);
    vCR(~isMatch) = 1;

    averaging_coefficients = computeJ1vCR(s4e,pos4n,n4sDb,counts_per_node,vCR);
    bubble_coefficients   = computeJ2vCR(n4s,vCR,averaging_coefficients);

    actSolution = integrateJ3vCR_lineload(...
        c4n,n4s_line,sides_line,averaging_coefficients,bubble_coefficients,lineload);

    verifyEqual(testCase, expSolution, actSolution, 'RelTol', 1e-6);
end

