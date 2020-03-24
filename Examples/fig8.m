% Figure 8: Dose-volume histogram for one PTV and one OAR with one 
% dose-volume constraint

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
rectum.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};

% Create problem instance
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);
xMat = prob.x0;
disp('Initialization')
fprintf('OAR %% > 50 Gy: %.2f, PTV D95: %.2f\n\n',prob.getPercent(2,1,xMat),...
    prob.getPercentile(prob.structs{1}.A*xMat,0.95));
% OAR % > 50 Gy: 56.80, PTV D95: 79.65

% Load approximate dose
load(['ex1Results/ex1aApprox.mat'])
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose')
fprintf('OAR %% > 50 Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',...
    prob.getPercent(2,1,x),prob.getPercentile(prob.structs{1}.A*x,0.95),t);
% OAR % > 50 Gy: 52.73, PTV D95: 79.67, Time: 2.48

% Load approximate dose with continuation
load(['ex1Results/ex1aContinue.mat'])
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose with continuation')
fprintf('OAR %% > 50 Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',...
    prob.getPercent(2,1,x),prob.getPercentile(prob.structs{1}.A*x,0.95),t);
% OAR % > 50 Gy: 50.85, PTV D95: 79.67, Time: 65.70

% Load polished dose
load(['ex1Results/ex1aPolish.mat'])
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Polished dose')
fprintf('OAR %% > 50 Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',...
    prob.getPercent(2,1,x),prob.getPercentile(prob.structs{1}.A*x,0.95),t);
% OAR % > 50 Gy: 49.76, PTV D95: 79.65, Time: 7.46

% Plot dose-volume histograms
prob.plotDVHPaper(xMat)
