function tests = test_CRPoisson_lineload
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

function testRHS(testCase)
    b = CRPoisson_lineload(4)
end
    