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

% Load approximate dose
load('ex2Results/ex2Approx.mat')
x0 = results.x0;
x1 = results.x;
disp('Approximate dose')
fprintf('OAR %% > 20 Gy: %.2f, %% > 50 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,x1),prob.getPercent(2,2,x1),prob.getPercent(2,3,x1));
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x1,0.95),...
    prob.getPercent(1,2,x1),prob.getPercent(1,3,x1));

% Load polished dose
load('ex2Results/ex2Polish.mat')
x2 = results.x;
disp('Polished Dose')
fprintf('OAR %% > 20 Gy: %.2f, %% > 50 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,x2),prob.getPercent(2,2,x2),prob.getPercent(2,3,x2));
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x2,0.95),...
    prob.getPercent(1,2,x2),prob.getPercent(1,3,x2));

% Plot dose-volume histograms
prob.plotDVHPaper([x0 x1 x2])
