% Figure 12: Dose-volume histograms for prostate tumor with uniform dose
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
    struct('type','ldvc','dose',81,'percent',5,'weight',100),...
    struct('type','udvc','dose',85,'percent',0,'weight',100)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',20,'percent',60,'weight',1),...
    struct('type','udvc','dose',40,'percent',40,'weight',1),...
    struct('type','udvc','dose',60,'percent',20,'weight',1)};

% Create problem instance
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);
x0 = prob.x0;
disp('Initialization')
fprintf('OAR %% > 20 Gy: %.2f, %% > 40 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,x0),prob.getPercent(2,2,x0),prob.getPercent(2,3,x0));
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x0,0.95),...
    prob.getPercent(1,2,x0),prob.getPercent(1,3,x0));
% OAR % > 20 Gy: 68.81, % > 40 Gy: 61.29, % > 60 Gy: 38.470.
% PTV D95: 79.65, % < 81 Gy: 53.26, % > 85 Gy: 0.00, Time: 0.1792

% Load approximate dose
load('ex2Results/ex2Approx.mat')
x1 = results.x;
t1 = results.time;
disp('Approximate dose')
fprintf('OAR %% > 20 Gy: %.2f, %% > 40 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,x1),prob.getPercent(2,2,x1),prob.getPercent(2,3,x1));
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f, Time: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x1,0.95),...
    prob.getPercent(1,2,x1),prob.getPercent(1,3,x1),t1);
% OAR % > 20 Gy: 61.10, % > 40 Gy: 42.66, % > 60 Gy: 21.84
% PTV D95: 80.72, % < 81 Gy: 12.53, % > 85 Gy: 0.00, Time: 18.87

% Plot dose-volume histograms
prob.plotDVHPaper([x0 x1 x1])
