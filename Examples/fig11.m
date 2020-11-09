% Figure 11: Dose-volume histograms for multiple PTVs and OARs.

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
disp('Initialization');
fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 30 Gy: %.2f\n',...
    prob.getPercent(3,1,50),prob.getPercent(4,1,30));
fprintf('Prostate D95: %.2f, Lymph Nodes D95: %.2f\n\n',...
    prob.getPercentile(1,0.95),prob.getPercentile(2,0.95));
% Rectum % > 50 Gy: 82.47, Bladder % > 30 Gy: 91.87
% Prostate D95: 76.96, Lymph Nodes D95: 58.53
% Time: 0.368328

% Load approximate dose
load('ex3Results/ex3Approx.mat')
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose');
fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 30 Gy: %.2f\n',...
    prob.getPercent(3,1,50,x),prob.getPercent(4,1,30,x));
fprintf('Prostate D95: %.2f, Lymph Nodes D95: %.2f\n',...
    prob.getPercentile(1,0.95,x),prob.getPercentile(2,0.95,x));
calcTime = results.time;
fprintf('Time: %.2f\n\n',t);
% Rectum % > 50 Gy: 57.33, Bladder % > 30 Gy: 38.14
% Prostate D95: 76.70, Lymph Nodes D95: 57.51
% Time: 15.15 (total 15.52)

% Load approximate dose with continuation
load('ex3Results/ex3Continue.mat')
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose with continuation')
fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 30 Gy: %.2f\n',...
    prob.getPercent(3,1,50,x),prob.getPercent(4,1,30,x));
fprintf('Prostate D95: %.2f, Lymph Nodes D95: %.2f\n',...
    prob.getPercentile(1,0.95,x),prob.getPercentile(2,0.95,x));
calcTime = results.time;
fprintf('Time: %.2f\n',t);
% Rectum % > 50 Gy: 47.53, Bladder % > 30 Gy: 29.93
% Prostate D95: 76.48, Lymph Nodes D95: 57.31
% Time: 84.46 (total 84.83)

% Plot dose-volume histograms
prob.plotDVHPaper(xMat)
