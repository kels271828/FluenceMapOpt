% Figure 13: Dose-volume histograms for polished solutions for problem
% with multiple PTVs and OARs.

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
structs = {prostate,nodes,rectum,bladder};
prob = FluenceMapOpt(structs);
xMat = [];

% Load approximate dose
labels = ['a' 'd' 'b' 'f' 'e' 'c'];
for ii = 1:length(labels)
    fprintf('\nMethod %s\n',labels(ii));
    load(['ex3Results/ex3' labels(ii) 'Polish.mat'])
    x = results.x;
    xMat = [xMat x];
    fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 30 Gy: %.2f\n',...
        prob.getPercent(3,1,x),prob.getPercent(4,1,x));
    fprintf('Prostate D95: %.2f, Lymph Nodes D95: %.2f\n',...
        prob.getPercentile(prob.structs{1}.A*x,0.95),...
        prob.getPercentile(prob.structs{2}.A*x,0.95));
    prob.x = x;
    fprintf('Uniform Objective: %.2f\n',prob.getObj('unif'));
    calcTime = results.time;
    fprintf('Time: %.2f\n',calcTime);
end

% Constraint Generation
% Rectum % > 50 Gy: 42.79, Bladder % > 30 Gy: 27.50
% Prostate D95: 76.59, Lymph Nodes D95: 55.73
% Uniform Objective: 13.31
% Time: 40.03

% Iterative Method
% Rectum % > 50 Gy: 49.12, Bladder % > 30 Gy: 29.19
% Prostate D95: 76.62, Lymph Nodes D95: 55.55
% Uniform Objective: 13.13
% Time: 43.86

% Our Method
% Rectum % > 50 Gy: 49.91, Bladder % > 30 Gy: 29.98
% Prostate D95: 76.61, Lymph Nodes D95: 57.07
% Uniform Objective: 11.90
% Time: 39.31

% Our Method with Continuation
% Rectum % > 50 Gy: 49.91, Bladder % > 30 Gy: 30.00
% Prostate D95: 76.58, Lymph Nodes D95: 57.30
% Uniform Objective: 11.81
% Time: 38.59

% Slack Method
% Rectum % > 50 Gy: 49.73, Bladder % > 30 Gy: 29.91
% Prostate D95: 75.82, Lymph Nodes D95: 57.61
% Uniform Objective: 12.26
% Time: 34.50

% Convex Method
% Rectum % > 50 Gy: 38.34, Bladder % > 30 Gy: 26.38
% Prostate D95: 75.21, Lymph Nodes D95: 57.49
% Uniform Objective: 13.12
% Time: 35.25

% Plot dose-volume histograms
prob.plotDVHPaper(xMat)
