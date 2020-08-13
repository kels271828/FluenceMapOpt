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

% Calculate approximate dose
fprintf('\nCalculating approximate dose\n\n');
prob.calcBeams();
fprintf('\nIterations: %d, Time: %.2f\n',prob.nIter,prob.time);
prob.saveResults('ex1Approx.mat');

% Calculate approximate dose with continuation
fprintf('\nCalculating approximate dose with continuation 1\n\n');
[prob,nIter,ttime] = calcBeamsContinue(prob,structs,0.99,0.01,0,100,1,0);
fprintf('\nIterations: %d, Time: %.2f\n',nIter,ttime);
prob.saveResults('ex1Continue1.mat');

% Calculate approximate dose with continuation
fprintf('\nCalculating approximate dose with continuation 2\n\n');
prob = FluenceMapOpt(structs);
[prob,nIter,ttime] = calcBeamsContinue(prob,structs,0.99,0.01,1,100,1,0);
fprintf('\nIterations: %d, Time: %.2f\n',nIter,ttime);
prob.saveResults('ex1Continue2.mat');
