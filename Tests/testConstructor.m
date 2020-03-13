function tests = testConstructor
    currentFolder = pwd;
    cd ..
    addpath(genpath(pwd));
    cd(currentFolder);
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    unif.name = 'PTV_68';
    unif.terms = {struct('type','unif','dose',80,'weight',1)};
    ldvc.name = 'PTV_56';
    ldvc.terms = {struct('type','ldvc','dose',60,'percent',15,'weight',1)};
    udvc.name = 'Rectum';
    udvc.terms = {struct('type','udvc','dose',30,'percent',15,'weight',1)};
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
    assert(prob.overlap == false,'Incorrect overlap')
    assert(prob.tol == 1e-3,'Incorrect tol')
    assert(prob.maxIter == 1000,'Incorrect maxIter')
    assert(prob.nStructs == 3,'Incorrect number of structures')
    assert(prob.nAngles == 7,'Incorrect number of angles')
    assert(strcmp(prob.nnls,'minConf_TMP'),'Incorrect nnls solver')
end

function testVararginParameters(testCase)
    structs = testCase.TestData.structs;
    prob = FluenceMapOpt(structs,'maxIter',0,'tol',1,'x0',2,'lambda',3,...
        'overlap',4,'angles',6);
    assert(prob.angles == 6,'Incorrect angles')
    assert(prob.overlap == 4,'Incorrect overlap')
    assert(prob.lambda == 3,'Incorrect lambda')
    assert(prob.x0 == 2,'Incorrect x0')
    assert(prob.tol == 1,'Incorrect tol')
    assert(prob.maxIter == 0,'Incorrect maxIter')
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
    [A,~,~,~] = prob.getA('full');
    assert(size(A,1) == 17882,'Incorrect rows')
    assert(size(A,2) == prob.nBeamlets,'Incorrect columns')
    [A,~,~,~] = prob.getA('unif');
    assert(size(A,1) == 7756,'Incorrect rows')
    assert(size(A,2) == prob.nBeamlets,'Incorrect columns')
    [A,~,~,~] = prob.getA('slack');
    assert(size(A,1) == 17882,'Incorrect rows')
    assert(size(A,2) == 11112,'Incorrect columns')
end

function testGetD(testCase)
    prob = testCase.TestData.prob;
    d = prob.getd('unif');
    assert(length(d) == 7756,'Incorrect length')
end

function testX0(testCase)
    prob = testCase.TestData.prob;
    x0 = prob.projX('unif');
    percentUnder = sum(prob.structs{2}.A*x0 < 60)/prob.structs{2}.nVoxels;
    assert(abs(percentUnder - 0.8782) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{3}.A*x0 > 30)/prob.structs{3}.nVoxels;
    assert(abs(percentOver - 0.6427) < 1e-4,'Incorrect initialization')
end
