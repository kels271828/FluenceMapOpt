% Figure 11: Dose-volume histograms for prostate tumor with uniform dose
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
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x0,0.95),...
    prob.getPercent(1,2,x0),prob.getPercent(1,3,x0));
% OAR % > 20 Gy: 68.81, % > 50 Gy: 56.80, % > 60 Gy: 38.47
% PTV D95: 79.65, % < 81 Gy: 53.26, % > 85 Gy: 0.00

% Load approximate dose
load('ex2Results/ex2Approx.mat')
x1 = results.x;
t1 = results.time;
disp('Approximate dose')
fprintf('OAR %% > 20 Gy: %.2f, %% > 50 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,x1),prob.getPercent(2,2,x1),prob.getPercent(2,3,x1));
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f, Time: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x1,0.95),...
    prob.getPercent(1,2,x1),prob.getPercent(1,3,x1),t1);
% OAR % > 20 Gy: 61.71, % > 50 Gy: 47.88, % > 60 Gy: 23.54
% PTV D95: 79.66, % < 81 Gy: 43.40, % > 85 Gy: 0.00, Time: 10.31

% Load approximate dose with continuation
load('ex2Results/ex2Continue.mat')
x2 = results.x;
t2 = results.time;
disp('Approximate dose with continuation')
fprintf('OAR %% > 20 Gy: %.2f, %% > 50 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,x2),prob.getPercent(2,2,x2),prob.getPercent(2,3,x2));
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f, Time: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x2,0.95),...
    prob.getPercent(1,2,x2),prob.getPercent(1,3,x2),t2);
% OAR % > 20 Gy: 60.38, % > 50 Gy: 43.51, % > 60 Gy: 20.63
% PTV D95: 80.81, % < 81 Gy: 8.08, % > 85 Gy: 0.00, Time: 702.85

% Load polished dose
load('ex2Results/ex2Polish.mat')
x3 = results.x;
t3 = results.time;
disp('Polished Dose')
fprintf('OAR %% > 20 Gy: %.2f, %% > 50 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,x3),prob.getPercent(2,2,x3),prob.getPercent(2,3,x3));
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f, Time: %.2f\n\n',...
    prob.getPercentile(prob.structs{1}.A*x3,0.95),...
    prob.getPercent(1,2,x3),prob.getPercent(1,3,x3),t3);
% OAR % > 20 Gy: 59.28, % > 50 Gy: 36.89, % > 60 Gy: 19.90
% PTV D95: 81.00, % < 81 Gy: 4.76, % > 85 Gy: 0.00, Time: 145.44

% Plot dose-volume histograms
prob.plotDVHPaper([x0 x1 x2 x3])
