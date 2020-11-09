% Example 1: One PTV and one OAR with one dose-volume constraint

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};

% Create problem instance
fprintf('\nExample 1\n');
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);
x0 = prob.x0;

% Calculate approximate dose
fprintf('\nCalculating approximate dose\n\n');
prob.calcBeams();
fprintf('\nIterations: %d, Time: %.2f\n',prob.nIter,prob.time);
prob.saveResults('ex1Approx.mat');

% Calculate approximate dose with continuation
fprintf('\nCalculating approximate dose with continuation (a)\n\n');
prob.updateStructs(structs,x0);
prob = calcBeamsContinue(prob,structs,0.99,0.01,100,true,0,false);
fprintf('Iterations: %d, Time: %.2f\n',prob.nIter,prob.time);
prob.saveResults('ex1ContinueA.mat');

% Calculate approximate dose with continuation
fprintf('\nCalculating approximate dose with continuation (b)\n\n');
prob.updateStructs(structs,x0);
prob = calcBeamsContinue(prob,structs,0.99,0.01,100,true,1,false);
fprintf('Iterations: %d, Time: %.2f\n',prob.nIter,prob.time);
prob.saveResults('ex1ContinueB.mat');
