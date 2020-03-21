% Figure 10: Dose-volume histograms for prostate tumor with uniform dose
% target of 81 Gy with various dose-volume constraints on the rectum.

% NOTE: Modified FluenceMapOpt.plotConstraints to dotted line for uniform
% target dose.

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1),...
    struct('type','ldvc','dose',81,'percent',5,'weight',1),...
    struct('type','udvc','dose',85,'percent',0,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',20,'percent',60,'weight',1),...
    struct('type','udvc','dose',50,'percent',50,'weight',1),...
    struct('type','udvc','dose',60,'percent',20,'weight',1)};

% Create problem instance
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);
x0 = prob.x0;
disp('Initialization')
fprintf('OAR %% > 20 Gy: %.2f, %% > 50 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,x0),prob.getPercent(2,2,x0),prob.getPercent(2,3,x0));
% OAR % > 20 Gy: 68.81, % > 50 Gy: 56.80, % > 60 Gy: 38.47
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x0,0.95),...
    prob.getPercent(1,2,x0),prob.getPercent(1,3,x0));
% PTV D95: 79.65, % < 81 Gy: 53.26, % > 85 Gy: 0.00

% Load approximate dose
load('ex2Results/ex2Approx.mat')
x1 = results.x;
t1 = results.time;
disp('Approximate dose')
fprintf('OAR %% > 20 Gy: %.2f, %% > 50 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,x1),prob.getPercent(2,2,x1),prob.getPercent(2,3,x1));
% OAR % > 20 Gy: 61.35, % > 50 Gy: 45.27, % > 60 Gy: 23.24
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f, Time: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x1,0.95),...
    prob.getPercent(1,2,x1),prob.getPercent(1,3,x1),t1);
% PTV D95: 79.90, % < 81 Gy: 55.89, % > 85 Gy: 0.00, Time: 103.41

% Load polished dose
load('ex2Results/ex2Polish.mat')
x2 = results.x;
t2 = results.time;
disp('Polished Dose')
fprintf('OAR %% > 20 Gy: %.2f, %% > 50 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,x2),prob.getPercent(2,2,x2),prob.getPercent(2,3,x2));
% OAR % > 20 Gy: 59.89, % > 50 Gy: 36.59, % > 60 Gy: 19.96
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f, Time: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x2,0.95),...
    prob.getPercent(1,2,x2),prob.getPercent(1,3,x2),t2);
% PTV D95: 81.00, % < 81 Gy: 4.80, % > 85 Gy: 0.00, Time: 132.17

% Plot dose-volume histograms
prob.plotDVHPaper([x0 x1 x2])
