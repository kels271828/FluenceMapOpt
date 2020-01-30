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
    testCase.TestData.structs = {unif,ldvc,udvc};
    testCase.TestData.prob = FluenceMapOpt(testCase.TestData.structs);
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
    assert(prob.nStructs == 3,'Incorrect number of structures')
    assert(prob.nAngles == 7,'Incorrect number of angles')
end

function testVararginParameters(testCase)
    structs = testCase.TestData.structs;
    prob = FluenceMapOpt(structs,[],0,1,2,3,4);
    assert(prob.angles == 0,'Incorrect angles')
    assert(prob.lambda == 1,'Incorrect lambda')
    assert(prob.maxIter == 2,'Incorrect maxIter')
    assert(prob.tol == 3,'Incorrect tol')
    assert(prob.overlap == 4,'Incorrect overlap')
end

function testSizeD(testCase)
    prob = testCase.TestData.prob;
    [D,nBeamlts] = FluenceMapOpt.getD(prob.angles);
    assert(size(D,1) == 3047040,'Incorrect number of rows')
    assert(size(D,2) == 986,'Incorrect number of columns')
    assert(nBeamlts == 986,'Incorrect number of beamlts')
end
 
function testGetStructVars(testCase)
    prob = testCase.TestData.prob;
    for ii = 1:prob.nStructs
        assert(prob.structs{ii}.nVoxels == size(prob.structs{ii}.A,1),...
            'Incorrect number of voxels')
        assert(prob.structs{ii}.nTerms == length(prob.structs{ii}.terms),...
            'Incorrect number of terms')
    end
end

function testGetTermVars(testCase)
    prob = testCase.TestData.prob;
    for ii = 2:3
        nVoxels = prob.structs{ii}.nVoxels;
        percent = prob.structs{ii}.terms{1}.percent;
        k = prob.structs{ii}.terms{1}.k;
        assert(k <= percent*nVoxels/100,'Incorrect k')
    end
end
 
function testGetNames(testCase)
    prob = testCase.TestData.prob;
    names = FluenceMapOpt.getNames(prob.structs,prob.nStructs);
    for ii = 1:prob.nStructs
        assert(strcmp(names{ii},prob.structs{ii}.name),'Incorrect name');
    end
end

function testSizeA(testCase)
    prob = testCase.TestData.prob;
    [A,~,~] = FluenceMapOpt.getA(prob.structs,prob.lambda,prob.nStructs,...
        prob.nBeamlts);
    assert(size(A,1) == 17882,'Incorrect rows')
    assert(size(A,2) == prob.nBeamlts,'Incorrect columns')
    [A,~,~] = FluenceMapOpt.getA(prob.structs,prob.lambda,prob.nStructs,...
        prob.nBeamlts,'full');
    assert(size(A,1) == 17882,'Incorrect rows')
    assert(size(A,2) == prob.nBeamlts,'Incorrect columns')
    [A,~,~] = FluenceMapOpt.getA(prob.structs,prob.lambda,prob.nStructs,...
        prob.nBeamlts,'unif');
    assert(size(A,1) == 7756,'Incorrect rows')
    assert(size(A,2) == prob.nBeamlts,'Incorrect columns')
end

function testGetD(testCase)
    prob = testCase.TestData.prob;
    d = FluenceMapOpt.getd(prob.structs,prob.lambda,prob.nStructs,...
        prob.nBeamlts,'unif');
    assert(length(d) == 7756,'Incorrect length')
end

function testInitX(testCase)
    prob = testCase.TestData.prob;
    [A,lb,ub] = FluenceMapOpt.getA(prob.structs,prob.lambda,prob.nStructs,...
        prob.nBeamlts,'unif');
    d = FluenceMapOpt.getd(prob.structs,prob.lambda,prob.nStructs,...
        prob.nBeamlts,'unif');
    xInit = FluenceMapOpt.getInitX(A,d,lb,ub);
    percentOver = sum(prob.structs{3}.A*xInit > 30)/prob.structs{3}.nVoxels;
    assert(abs(percentOver - 0.6470) < 1e-2)
end
