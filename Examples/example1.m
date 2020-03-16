% One PTV and one OAR with one dose-volume constraint

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% Loop over different dose/percent values
labels = ['a' 'b' 'c'];
doses = [50 30 10];
for ii = 1:length(doses)
    fprintf('\nExample 1%s\n',labels(ii));
    
    % PTV - prostate
    prostate.name = 'PTV_68';
    prostate.terms = {struct('type','unif','dose',81,'weight',1)};

    % OAR - rectum
    rectum.name = 'Rectum';
    rectum.terms = {struct('type','udvc','dose',doses(ii),...
        'percent',doses(ii),'weight',1)};

    % Create problem instance
    structs = {prostate,rectum};
    prob = FluenceMapOpt(structs,'tol',5e-5);

    % Calculate approximate dose
    fprintf('\nCalculating approximate dose\n\n');
    prob.calcBeams();
    fprintf('\nIterations: %d, Time: %.2f\n',prob.nIter,prob.time);
    filename1 = ['e1' labels(ii) '_x1.mat'];
    prob.saveResults(filename1);

    % Calculate polished dose
    fprintf('\nCalculating polished dose\n');
    prob.calcBeamsPolish(prob.x);
    fprintf('\nTime: %.2f\n',prob.time);
    filename2 = ['e1' labels(ii) '_x2.mat'];
    prob.saveResults(filename2);

    % Calculate approximate dose with continuation
    fprintf('\nCalculating approximate dose with continuation\n\n');
    prob = FluenceMapOpt(structs,'tol',1e-3);
    prob = calcBeamsCont(prob,structs,0.5,1.5,1e-3,100,true);
    fprintf('\nIterations: %d, Time: %.2f\n',prob.nIter,prob.time);
    filename3 = ['e1' labels(ii) '_x3.mat'];
    prob.saveResults(filename3);
end
