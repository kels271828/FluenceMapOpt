% Figure 12: Dose-volume histograms for approximate solutions for problem
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
xMat = prob.x0;
fprintf('\nMethod a\n');
fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 30 Gy: %.2f\n',...
    prob.getPercent(3,1,prob.x0),prob.getPercent(4,1,prob.x0));
fprintf('Prostate D95: %.2f, Lymph Nodes D95: %.2f\n',...
    prob.getPercentile(prob.structs{1}.A*prob.x0,0.95),...
    prob.getPercentile(prob.structs{2}.A*prob.x0,0.95));
prob.x = prob.x0;
fprintf('Uniform Objective: %.2f\n',prob.getObj('unif'));

% Load approximate dose
labels = ['d' 'b' 'f' 'e' 'c'];
for ii = 1:length(labels)
    fprintf('\nMethod %s\n',labels(ii));
    load(['ex3Results/ex3' labels(ii) 'Approx.mat'])
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
% Rectum % > 50 Gy: 82.47, Bladder % > 30 Gy: 91.87
% Prostate D95: 76.96, Lymph Nodes D95: 58.53
% Uniform Objective: 10.49
% Time: 0.323528

% Iterative Method
% Rectum % > 50 Gy: 65.67, Bladder % > 30 Gy: 58.14
% Prostate D95: 75.70, Lymph Nodes D95: 55.23
% Uniform Objective: 12.38
% Time: 3.37

% Our Method
% Rectum % > 50 Gy: 57.33, Bladder % > 30 Gy: 38.14
% Prostate D95: 76.70, Lymph Nodes D95: 57.51
% Uniform Objective: 11.39
% Time: 15.42

% Our Method with Continuation
% Rectum % > 50 Gy: 51.67, Bladder % > 30 Gy: 30.82
% Prostate D95: 76.54, Lymph Nodes D95: 57.27
% Uniform Objective: 11.82
% Time: 1231.80

% Slack Method
% Rectum % > 50 Gy: 53.38, Bladder % > 30 Gy: 35.48
% Prostate D95: 75.89, Lymph Nodes D95: 57.75
% Uniform Objective: 12.09
% Time: 1370.15

% Convex Method
% Rectum % > 50 Gy: 19.84, Bladder % > 30 Gy: 12.37
% Prostate D95: 71.98, Lymph Nodes D95: 54.52
% Uniform Objective: 23.45
% Time: 267.78

% Plot dose-volume histograms
prob.plotDVHPaper(xMat)
