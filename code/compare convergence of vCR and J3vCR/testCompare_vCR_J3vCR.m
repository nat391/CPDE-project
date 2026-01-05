

%% Test 1: testing right-hand side f = 0
[x, x_J3] = compare_vCR_J3vCR(4,@(x) 0);
assert(isequal(x, zeros(size(x))))
assert(isequal(x_J3, zeros(size(x_J3))))
