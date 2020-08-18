% Figure 10: Dose-volume histograms for prostate tumor with uniform dose
% target of 81 Gy with various dose-volume constraints on the rectum.

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1),...
    struct('type','ldvc','dose',81,'percent',5,'weight',10),...
    struct('type','udvc','dose',85,'percent',0,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',20,'percent',60,'weight',1),...
    struct('type','udvc','dose',40,'percent',40,'weight',1),...
    struct('type','udvc','dose',60,'percent',20,'weight',1)};

% Create problem instance
fprintf('Example2\n\n')
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);
xMat = prob.x0;
disp('Initialization')
fprintf('OAR %% > 20 Gy: %.2f, %% > 40 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,20),prob.getPercent(2,2,40),...
    prob.getPercent(2,3,60));
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f\n\n',...
    prob.getPercentile(1,0.95,xMat),...
    prob.getPercent(1,2,81,xMat),prob.getPercent(1,3,85,xMat));
% OAR % > 20 Gy: 68.81, % > 40 Gy: 61.29, % > 60 Gy: 38.47
% PTV D95: 79.65, % < 81 Gy: 53.26, % > 85 Gy: 0.00, Time: 0.1792

% Load approximate dose
load('ex2Results/ex2Approx.mat')
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose')
fprintf('OAR %% > 20 Gy: %.2f, %% > 40 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,20,x),prob.getPercent(2,2,40,x),...
    prob.getPercent(2,3,60,x));
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f, Time: %.2f\n\n',...
    prob.getPercentile(1,0.95,x),...
    prob.getPercent(1,2,81,x),prob.getPercent(1,3,85,x),t);
% OAR % > 20 Gy: 61.95, % > 40 Gy: 42.90, % > 60 Gy: 21.84
% PTV D95: 80.26, % < 81 Gy: 24.27, % > 85 Gy: 0.00, Time: 8.71 (total 8.89)

% Load approximate dose with continuation
load('ex2Results/ex2Continue.mat')
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose with continuation')
fprintf('OAR %% > 20 Gy: %.2f, %% > 40 Gy: %.2f, %% > 60 Gy: %.2f\n',...
    prob.getPercent(2,1,20,x),prob.getPercent(2,2,40,x),...
    prob.getPercent(2,3,60,x));
fprintf('PTV D95: %.2f, %% < 81 Gy: %.2f, %% > 85 Gy: %.2f, Time: %.2f\n\n',...
    prob.getPercentile(1,0.95,x),...
    prob.getPercent(1,2,81,x),prob.getPercent(1,3,85,x),t);
% OAR % > 20 Gy: 59.10, % > 40 Gy: 39.26, % > 60 Gy: 19.78
% PTV D95: 81.32, % < 81 Gy: 4.09, % > 85 Gy: 0.00, Time: 55.69 (total 55.86)

% Plot dose-volume histograms
prob.plotDVHPaper(xMat,false)
