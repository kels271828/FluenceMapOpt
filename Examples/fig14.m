% Figure 14: Dose-volume histograms for multiple PTVs and OARs.

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
fprintf('Example 3\n\n')
structs = {prostate,nodes,rectum,bladder};
prob = FluenceMapOpt(structs);
xMat = prob.x0;
fprintf('\nInitialization\n');
fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 30 Gy: %.2f\n',...
    prob.getPercent(3,1,prob.x0),prob.getPercent(4,1,prob.x0));
fprintf('Prostate D95: %.2f, Lymph Nodes D95: %.2f\n',...
    prob.getPercentile(prob.structs{1}.A*prob.x0,0.95),...
    prob.getPercentile(prob.structs{2}.A*prob.x0,0.95));
% Rectum % > 50 Gy: 82.47, Bladder % > 30 Gy: 91.87
% Prostate D95: 76.96, Lymph Nodes D95: 58.53
% Time: 0.368328

% Load approximate dose
fprintf('\nApproximate dose\n');
load(['ex3Results/ex3Approx.mat'])
x = results.x;
xMat = [xMat x];
fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 30 Gy: %.2f\n',...
    prob.getPercent(3,1,x),prob.getPercent(4,1,x));
fprintf('Prostate D95: %.2f, Lymph Nodes D95: %.2f\n',...
    prob.getPercentile(prob.structs{1}.A*x,0.95),...
    prob.getPercentile(prob.structs{2}.A*x,0.95));
calcTime = results.time;
fprintf('Time: %.2f\n',calcTime);
% Rectum % > 50 Gy: 57.33, Bladder % > 30 Gy: 38.14
% Prostate D95: 76.70, Lymph Nodes D95: 57.51
% Time: 14.74

% Plot dose-volume histograms
prob.plotDVHPaper(xMat)
