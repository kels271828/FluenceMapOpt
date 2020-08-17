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
    struct('type','ldvc','dose',81,'percent',5,'weight',10),...
    struct('type','udvc','dose',85,'percent',0,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',20,'percent',60,'weight',1),...
    struct('type','udvc','dose',40,'percent',40,'weight',1),...
    struct('type','udvc','dose',60,'percent',20,'weight',1)};

% Create problem instance
fprintf('\nExample 2\n');
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);
x0 = prob.x0;

% Calculate approximate dose
fprintf('\nCalculating approximate dose\n\n');
prob = calcBeamsConsolidate(prob,true);
fprintf('Iterations: %d, Time: %.2f\n',prob.nIter,prob.time);
prob.saveResults('ex2Approx.mat');

% Calculate approximate dose with continuation
fprintf('\nCalculating approximate dose with continuation\n\n');
prob.updateStructs(structs,x0);
prob = calcBeamsContinue(prob,structs,0.99,0.01,100,true,0,true);
fprintf('Iterations: %d, Time: %.2f\n',prob.nIter,prob.time);
prob.saveResults('ex2Continue.mat');
