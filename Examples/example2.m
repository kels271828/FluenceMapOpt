% Example 2: One PTV and one OAR with multiple dose-volume constraints

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1),...
    struct('type','ldvc','dose',81,'percent',5,'weight',1),...
    struct('type','udvc','dose',85,'percent',0,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',20,'percent',60,'weight',1),...
    struct('type','udvc','dose',50,'percent',50,'weight',1),...
    struct('type','udvc','dose',60,'percent',20,'weight',1)};

% Create problem instance
fprintf('Example 2\n');
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);

% Calculate approximate dose
fprintf('\nCalculating approximate dose\n\n');
prob.calcBeams();
fprintf('\nIterations: %d, Time: %.2f\n',prob.nIter,prob.time);
prob.saveResults('ex2Approx.mat');

% Calculate polished dose
fprintf('\nCalculating polished dose\n');
prob.calcBeamsPolish(prob.x);
fprintf('\nTime: %.2f\n',prob.time);
prob.saveResults('ex2Polish.mat');
