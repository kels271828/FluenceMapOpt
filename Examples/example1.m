% Example 1: One PTV and one OAR with one dose-volume constraint

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
    filename1 = ['ex1' labels(ii) 'Approx.mat'];
    prob.saveResults(filename1);
end
