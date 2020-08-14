% Figure 8: Convergence of objective function for one PTV and one OAR with
% one dose-volume constraint.

% Note: Change SetAccess of FluenceMapOpt properties structs, obj, wDiff,
% and nIter.

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
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);
prob.loadResults('ex1Results/ex1Approx.mat');

% Plot convergence
prob.plotObjPaper();
