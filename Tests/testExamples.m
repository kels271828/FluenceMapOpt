function tests = testExamples
    currentFolder = pwd;
    cd ..
    addpath(genpath(pwd));
    cd(currentFolder);
    tests = functiontests(localfunctions);
end

function testExample1(~)
    unif.name = 'PTV_68';
    unif.terms = {struct('type','unif','dose',81,'weight',1)};
    udvc.name = 'Rectum';
    udvc.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};
    prob = FluenceMapOpt({unif,udvc});
    percentOver = sum(prob.structs{2}.A*prob.x0 > 50)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.5680) < 1e-4,'Incorrect initialization')
    prob.calcBeams(false);
    assert(prob.nIter == 11,'Wrong number of iterations')
    percentOver = sum(prob.structs{2}.A*prob.x > 50)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.5273) < 1e-4,'Incorrect solution')
end

function testExample2(~)
    unif.name = 'PTV_68';
    unif.terms = {struct('type','unif','dose',81,'weight',1)};
    udvc.name = 'Rectum';
    udvc.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};
    prob = FluenceMapOpt({unif,udvc});
    percentOver = sum(prob.structs{2}.A*prob.x0 > 30)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.6414) < 1e-4,'Incorrect initialization')
    prob.calcBeams(false);
    assert(prob.nIter == 42,'Wrong number of iterations')
    percentOver = sum(prob.structs{2}.A*prob.x > 30)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.3416) < 1e-4,'Incorrect solution')
end

function testExample3(~)
    unif.name = 'PTV_68';
    unif.terms = {struct('type','unif','dose',81,'weight',1)};
    udvc.name = 'Rectum';
    udvc.terms = {struct('type','udvc','dose',10,'percent',10,'weight',1)};
    prob = FluenceMapOpt({unif,udvc});
    percentOver = sum(prob.structs{2}.A*prob.x0 > 10)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.7397) < 1e-4,'Incorrect initialization')
    prob.calcBeams(false);
    assert(prob.nIter == 30,'Wrong number of iterations')
    percentOver = sum(prob.structs{2}.A*prob.x > 10)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.2476) < 1e-4,'Incorrect solution')
end

function testExample4(~)
    unif.name = 'PTV_68';
    unif.terms = {struct('type','unif','dose',81,'weight',1),...
        struct('type','ldvc','dose',81','percent',5,'weight',1),...
        struct('type','udvc','dose',85,'percent',0,'weight',1)};
    udvc.name = 'Rectum';
    udvc.terms = {struct('type','udvc','dose',20,'percent',60,'weight',1),...
        struct('type','udvc','dose',40,'percent',40,'weight',1),...
        struct('type','udvc','dose',60,'percent',20,'weight',1)};

    % Check initialization
    prob = FluenceMapOpt({unif,udvc});
    percentUnder = sum(prob.structs{1}.A*prob.x0 < 81)/prob.structs{1}.nVoxels;
    assert(abs(percentUnder - 0.5326) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{1}.A*prob.x0 > 85)/prob.structs{1}.nVoxels;
    assert(abs(percentOver - 0) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{2}.A*prob.x0 > 20)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.6881) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{2}.A*prob.x0 > 40)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.6129) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{2}.A*prob.x0 > 60)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.3847) < 1e-4,'Incorrect initialization') 
   
    % Check solution
    prob.calcBeams(false);
    assert(prob.nIter == 32,'Wrong number of iterations')
    percentUnder = sum(prob.structs{1}.A*prob.x < 81)/prob.structs{1}.nVoxels;
    assert(abs(percentUnder - 0.5594) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{1}.A*prob.x > 85)/prob.structs{1}.nVoxels;
    assert(abs(percentOver - 0) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{2}.A*prob.x > 20)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.6280) < 1e-4,'Incorrect solution')
    percentOver = sum(prob.structs{2}.A*prob.x > 40)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.4399) < 1e-4,'Incorrect solution')
    percentOver = sum(prob.structs{2}.A*prob.x > 60)/prob.structs{2}.nVoxels;
    assert(abs(percentOver - 0.2330) < 1e-4,'Incorrect solution')
end

function testExample5(~)
    unif1.name = 'PTV_68';
    unif1.terms = {struct('type','unif','dose',81,'weight',1)};
    unif2.name = 'PTV_56';
    unif2.terms = {struct('type','unif','dose',60,'weight',1)};
    udvc1.name = 'Rectum';
    udvc1.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};
    udvc2.name = 'Bladder';
    udvc2.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};
    
    % Check initialization
    prob = FluenceMapOpt({unif1,unif2,udvc1,udvc2});
    percentOver = sum(prob.structs{3}.A*prob.x0 > 50)/prob.structs{3}.nVoxels;
    assert(abs(percentOver - 0.8247) < 1e-4,'Incorrect initialization')
    percentOver = sum(prob.structs{4}.A*prob.x0 > 30)/prob.structs{4}.nVoxels;
    assert(abs(percentOver - 0.9187) < 1e-4,'Incorrect initialization')
    
    % Check solution
    prob.calcBeams(false);
    assert(prob.nIter == 48,'Wrong number of iterations')
    percentOver = sum(prob.structs{3}.A*prob.x > 50)/prob.structs{3}.nVoxels;
    assert(abs(percentOver - 0.5733) < 1e-4,'Incorrect solution')
    percentOver = sum(prob.structs{4}.A*prob.x > 30)/prob.structs{4}.nVoxels;
    assert(abs(percentOver - 0.3814) < 1e-4,'Incorrect solution')
end