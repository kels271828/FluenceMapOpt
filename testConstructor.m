function tests = testConstructor
    addpath(genpath('PROSTATE'));
    addpath(genpath('Prostate_Dicom'));
    addpath(genpath('minConf'));
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    unif = struct('name','PTV_68');
    unif.terms = {struct('type','unif','dose',81,'weight',1)};
    ldvc = struct('name','PTV_56');
    ldvc.terms = {struct('type','ldvc','dose',60,'percent',60,'weight',1)};
    udvc = struct('name','Rectum');
    udvc.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};
    structs = {unif,ldvc,udvc};
    testCase.TestData.prob = FluenceMapOpt(structs);
end
    
function testDefaultParameters(testCase)
    prob = testCase.TestData.prob;
    defaultAngles = 0:52:328;
    for ii = 1:length(prob.angles)
        assert(prob.angles(ii) == defaultAngles(ii),'Incorrect angle')
    end
    assert(prob.lambda == 1e-8,'Incorrect lambda')
    assert(prob.maxIter == 500,'Incorrect maxIter')
    assert(prob.tol == 1e-3,'Incorrect tol')
    assert(prob.overlap == false,'Incorrect overlap')
    assert(size(prob.D,1) == 3047040,'Incorrect rows')
    assert(size(prob.D,2) == 986,'Incorrect columns')
end

function testVararginParameters(testCase)
    prob = testCase.TestData.prob;
    prob.setInputVars({[],0,1,2,3,4},6);
    assert(prob.angles == 0,'Incorrect angles')
    assert(prob.lambda == 1,'Incorrect lambda')
    assert(prob.maxIter == 2,'Incorrect maxIter')
    assert(prob.tol == 3,'Incorrect tol')
    assert(prob.overlap == 4,'Incorrect overlap')
end

function testSizeD(testCase)
    prob = testCase.TestData.prob;
    assert(size(prob.D,1) == 3047040,'Incorrect rows')
    assert(size(prob.D,2) == 986,'Incorrect columns')
end

function testGetTermVars(testCase)
    prob = testCase.TestData.prob;
    for ii = 2:3
        nVoxels = size(prob.structs{ii}.A,1);
        percent = prob.structs{ii}.terms{1}.percent;
        k = prob.structs{ii}.terms{1}.k;
        assert(k <= percent*nVoxels/100,'Incorrect k')
    end
end

function testSizeA(testCase)
    prob = testCase.TestData.prob;
    assert(size(prob.A,1) == 17882,'Incorrect rows')
    assert(size(prob.A,2) == 986,'Incorrect columns')
end

function testInitX(testCase)
    prob = testCase.TestData.prob;
    A = prob.structs{3}.A;
    nVoxels = size(A,1);
    percentOver = sum(A*prob.xInit > 30)/nVoxels;
    assert(abs(percentOver - 0.6470) < 1e-2)
end
