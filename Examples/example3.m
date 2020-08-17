% Example 3: Multiple PTVs and OARs

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% PTV - lymph nodes
nodes.name = 'PTV_56';
nodes.terms = {struct('type','unif','dose',60,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};

% OAR - bladder
bladder.name = 'Bladder';
bladder.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};

% Create problem instance
fprintf('\nExample3\n')
structs = {prostate,nodes,rectum,bladder};
prob = FluenceMapOpt(structs);
x0 = prob.x0;

% Calculate approximate dose
fprintf('\nCalculating approximate dose\n\n');
prob.calcBeams();
fprintf('\nIterations: %d, Time: %.2f\n\n',prob.nIter,prob.time);
prob.saveResults('ex3Approx.mat');

% Calculate approximate dose with continuation
fprintf('Calculating approximate dose with continuation\n\n');
prob.updateStructs(structs,x0);
prob = calcBeamsContinue(prob,structs,0.99,0.01,100,true,0,false);
fprintf('Iterations: %d, Time: %.2f\n',prob.nIter,prob.time);
prob.saveResults('ex3Continue.mat');
