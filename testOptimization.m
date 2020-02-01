function tests = testOptimization
    addpath(genpath('PROSTATE'));
    addpath(genpath('minConf'));
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    unif = struct('name','PTV_68');
    unif.terms = {struct('type','unif','dose',80,'weight',1)};
    ldvc = struct('name','PTV_56');
    ldvc.terms = {struct('type','ldvc','dose',60,'percent',15,'weight',1)};
    udvc = struct('name','Rectum');
    udvc.terms = {struct('type','udvc','dose',30,'percent',15,'weight',1)};
    testCase.TestData.structs = {unif,ldvc,udvc};
    testCase.TestData.prob = FluenceMapOpt(testCase.TestData.structs);
end

function testProjW(~)
    n = 100;
    k = 25;
    w1 = randn(n,1);
    w2 = FluenceMapOpt.projW(w1,k);
    assert(sum(w1 < 0) == sum(w2 < 0),'Wrong number negative')
    assert(sum(w2 <= 0) >= (n - k),'Wrong number nonpositive')
    assert(sum(w2 > 0) <= k,'Wrong number positive')
end

function testInitW(testCase)
    prob = testCase.TestData.prob;
    prob.initProb(false);
    for ii = 2:3
        k = prob.structs{ii}.terms{1}.k;
        w = prob.structs{ii}.terms{1}.w;
        assert(sum(w > 0) <= k,'Wrong projection') 
    end
end

function testInitObj(testCase)
    prob = testCase.TestData.prob;
    prob.initProb(false);
    assert(abs(prob.obj(1) - 2.2271e+03) < 1e0,'Wrong objective value')
end

function testUpdateX(testCase)
    prob = testCase.TestData.prob;
    prob.initProb(false);
    x = prob.projX('full');
    percentUnder = sum(prob.structs{2}.A*x < 60)/prob.structs{2}.nVoxels;
    assert(abs(percentUnder - 0.4828) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{3}.A*x > 30)/prob.structs{3}.nVoxels;
    assert(abs(percentOver - 0.4303) < 1e-4,'Incorrect initialization')
end

% test update w?
% test calcBeamlets()?
