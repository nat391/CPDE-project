function tests = test_integrateJ3vCR_slit
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    refinement_lvl = 2;
    [c4n,n4e,n4sDb] = loadGeometry('BigSquare',refinement_lvl);
    
    % find nodes that are on the slit
    nodes_slit = find(c4n(:,1)>=0 & c4n(:,2)==0);

    n4s = computeN4s(n4e);

    % get the indeces of the sides that are on the slit
    sides_slit = find( all( ismember(n4s, nodes_slit), 2 ) );

    % store everything
    testCase.TestData.c4n        = c4n;
    testCase.TestData.n4e        = n4e;
    testCase.TestData.n4s       = n4s;
    testCase.TestData.n4sDb     = n4sDb;
    testCase.TestData.sides_slit = sides_slit;
end


function testZeroInput(testCase)

    % retrieve shared data
    c4n        = testCase.TestData.c4n;
    n4e        = testCase.TestData.n4e;
    n4s        = testCase.TestData.n4s;
    n4sDb      = testCase.TestData.n4sDb;
    sides_slit = testCase.TestData.sides_slit;

    expSolution = 0;
    
    nr_sides = size(n4s,1);
    vCR = zeros(nr_sides,1);

    averaging_coefficients = computeJ1vCR(c4n,n4e,n4sDb,vCR);
    bubble_coefficients   = computeJ2vCR(n4s,vCR,averaging_coefficients);

    actSolution = integrateJ3vCR_slit( ...
        c4n,n4s,averaging_coefficients,bubble_coefficients,sides_slit);

    verifyEqual(testCase, expSolution, actSolution)
end


function testOneInput(testCase)

    % retrieve shared data
    c4n        = testCase.TestData.c4n;
    n4e        = testCase.TestData.n4e;
    n4s        = testCase.TestData.n4s;
    n4sDb      = testCase.TestData.n4sDb;
    sides_slit = testCase.TestData.sides_slit;
    
    nr_sides = size(n4s,1);
    
    % find dofs -> sides that are not on the boundary
    n4s = sort(n4s,2);
    n4sDb = sort(n4sDb,2);
    [isMatch, ~] = ismember(n4s,n4sDb,'rows');
    
    % define CR-input
    vCR = zeros(nr_sides);
    vCR(~isMatch) = 1;

    averaging_coefficients = computeJ1vCR(c4n,n4e,n4sDb,vCR);
    bubble_coefficients   = computeJ2vCR(n4s,vCR,averaging_coefficients);

    result = integrateJ3vCR_slit( ...
        c4n,n4s,averaging_coefficients,bubble_coefficients,sides_slit);

    % each node has 6 neighbouring triangles
    % for a fine mesh, most nodes along the slit will have 4 neighbouring
    % triangles with value 1 at that node. Therefore the nodal averaging
    % coefficient is 2/3. The bubble coefficient is 
    % vCR - J1vCR, so 1 - 2/3 = 1/3. The mean of both is 0.5.
    % Therefore, considering the boundary value, one expects a result that
    % tends to 0.5 with finer mesh size
    verifyThat(testCase, result, matlab.unittest.constraints.IsLessThan(0.5));
end

