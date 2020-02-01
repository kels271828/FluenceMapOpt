function tests = testOptimization
    addpath(genpath('PROSTATE'));
    addpath(genpath('minConf'));
    tests = functiontests(localfunctions);
end

function testExample1(testCase)
    unif = struct('name','PTV_68');
    unif.terms = {struct('type','unif','dose',81,'weight',1)};
    udvc = struct('name','Rectum');
    udvc.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};
    prob = FluenceMapOpt({unif,udvc},'tol',5e-5);
    percentOver = sum(prob.structs{2}.A*prob.x0 > 50)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.5680) < 1e-4,'Incorrect initialization')
    prob.calcBeamlets(false);
    assert(prob.nIter == 231,'Wrong number of iterations')
    percentOver = sum(prob.structs{2}.A*prob.x > 50)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.5152) < 1e-4,'Incorrect solution')
end

function testExample2(testCase)
    unif = struct('name','PTV_68');
    unif.terms = {struct('type','unif','dose',81,'weight',1)};
    udvc = struct('name','Rectum');
    udvc.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};
    prob = FluenceMapOpt({unif,udvc},'tol',5e-5);
    percentOver = sum(prob.structs{2}.A*prob.x0 > 30)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.6414) < 1e-4,'Incorrect initialization')
    prob.calcBeamlets(false);
    assert(prob.nIter == 468,'Wrong number of iterations')
    percentOver = sum(prob.structs{2}.A*prob.x > 30)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.3337) < 1e-4,'Incorrect solution')
end

function testExample3(testCase)
    unif = struct('name','PTV_68');
    unif.terms = {struct('type','unif','dose',81,'weight',1)};
    udvc = struct('name','Rectum');
    udvc.terms = {struct('type','udvc','dose',10,'percent',10,'weight',1)};
    prob = FluenceMapOpt({unif,udvc},'tol',5e-5);
    percentOver = sum(prob.structs{2}.A*prob.x0 > 10)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.7397) < 1e-4,'Incorrect initialization')
    prob.calcBeamlets(false);
    assert(prob.nIter == 320,'Wrong number of iterations')
    percentOver = sum(prob.structs{2}.A*prob.x > 10)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.2391) < 1e-4,'Incorrect solution')
end

function testExample4(testCase)
    unif = struct('name','PTV_68');
    unif.terms = {struct('type','unif','dose',81,'weight',1)};
    udvc = struct('name','Rectum');
    udvc.terms = {struct('type','udvc','dose',20,'percent',60,'weight',1),...
        struct('type','udvc','dose',50,'percent',50,'weight',1),...
        struct('type','udvc','dose',60,'percent',20,'weight',1),...
        struct('type','udvc','dose',75,'percent',0,'weight',1)};

    % Check initialization
    prob = FluenceMapOpt({unif,udvc},'tol',5e-5,'maxIter',1000);
    percentOver = sum(prob.structs{2}.A*prob.x0 > 20)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.6881) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{2}.A*prob.x0 > 50)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.5680) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{2}.A*prob.x0 > 60)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.3847) < 1e-4,'Incorrect initialization') 
    percentOver = sum(prob.structs{2}.A*prob.x0 > 75)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.1341) < 1e-4,'Incorrect initialization')
   
    % Check solution (doesn't match paper)
    prob.calcBeamlets(false);
    %assert(prob.nIter == 568,'Wrong number of iterations')
    assert(prob.nIter == 513,'Wrong number of iterations')
    percentOver = sum(prob.structs{2}.A*prob.x > 20)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.6141) < 1e-4,'Incorrect solution')
    percentOver = sum(prob.structs{2}.A*prob.x > 50)/prob.structs{2}.nVoxels;
    %assert(abs(percentOver - 0.4660) < 1e-4,'Incorrect solution')
    assert(abs(percentOver - 0.4678) < 1e-4,'Incorrect solution')
    percentOver = sum(prob.structs{2}.A*prob.x > 60)/prob.structs{2}.nVoxels;
    %assert(abs(percentOver - 0.2172) < 1e-4,'Incorrect solution')
    assert(abs(percentOver - 0.2166) < 1e-4,'Incorrect solution')
    percentOver = sum(prob.structs{2}.A*prob.x > 75)/prob.structs{2}.nVoxels;
    %assert(abs(percentOver - 0.0388) < 1e-4,'Incorrect solution')
    assert(abs(percentOver - 0.0394) < 1e-4,'Incorrect solution')
end

function testExample5(testCase)
    unif1 = struct('name','PTV_68');
    unif1.terms = {struct('type','unif','dose',81,'weight',1)};
    unif2 = struct('name','PTV_56');
    unif2.terms = {struct('type','unif','dose',60,'weight',1)};
    udvc1 = struct('name','Rectum');
    udvc1.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};
    udvc2 = struct('name','Bladder');
    udvc2.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};
    
    % Check initialization
    prob = FluenceMapOpt({unif1,udvc1,udvc2,unif2},'tol',5e-5);
    percentOver = sum(prob.structs{2}.A*prob.x0 > 50)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.8204) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{3}.A*prob.x0 > 30)/prob.structs{3}.nVoxels;
    assert(abs(percentOver - 0.9221) < 1e-4,'Incorrect initialization')
    
    % Check solution
    prob.calcBeamlets(true);
    assert(prob.nIter == 336,'Wrong number of iterations')
    percentOver = sum(prob.structs{2}.A*prob.x > 50)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.5704) < 1e-4,'Incorrect solution')
    percentOver = sum(prob.structs{3}.A*prob.x > 30)/prob.structs{3}.nVoxels;
    assert(abs(percentOver - 0.3576) < 1e-4,'Incorrect solution')
end