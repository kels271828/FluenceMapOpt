function tests = testOptimization
    addpath(genpath('PROSTATE'));
    addpath(genpath('Prostate_Dicom'));
    addpath(genpath('minConf'));
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    unif = struct('name','PTV_68');
    unif.terms = {struct('type','unif','dose',81,'weight',1)};
    ldvc = struct('name','PTV_56');
    udvc.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};
    structs = {unif,ldvc,udvc};
    testCase.TestData.prob = FluenceMapOpt(structs);
end

% Test project w (with udvc and ldvc)
% Test initW (with udvc and ldvc)

% align everything with our shiny new precise definitions