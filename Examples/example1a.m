% One PTV and one OAR with one dose-volume constraint

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% Set up problem

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};

% Create problem instance
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);

%% Calculate approximate dose

prob.calcBeams();
fprintf('\nIterations: %d, Time: %.2f\n\n',prob.nIter,prob.time);
prob.saveResults('e1a_x1.mat');

%% Calculate polished dose

prob.calcBeamsPolish(prob.x);
fprintf('\nTime: %.2f\n\n',prob.time);
prob.saveResults('e1a_x2.mat');

%% Calculate approximate dose with continuation

prob = calcBeamsCont(prob,structs,true);
fprintf('\nIterations: %d, Time: %.2f\n\n',prob.nIter,prob.time);
prob.saveResults('e1a_x3.mat');
