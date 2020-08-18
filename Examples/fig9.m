% Figures 9: Dose-volume histogram for one PTV and one OAR with one 
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
rectum.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};

% Create problem instance
fprintf('Example 1\n\n');
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);    
xMat = prob.x0;
disp('Initialization')
fprintf('OAR %% > 30 Gy: %.2f, PTV D95: %.2f\n\n',...
    prob.getPercent(2,1,30),prob.getPercentile(1,0.95));
% OAR % > 30 Gy: 64.14, PTV D95: 79.65, Time: 0.1792
    
% Load approximate dose
load('ex1Results/ex1Approx.mat')
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose')
fprintf('OAR %% > 30 Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',...
    prob.getPercent(2,1,30,x),prob.getPercentile(1,0.95,x),t);
% OAR % > 30 Gy: 34.16, PTV D95: 79.17, Time: 7.04 (total 7.22)

% Load approximate dose with continuation (a)
load('ex1Results/ex1ContinueA.mat')
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose with continuation (a)')
fprintf('OAR %% > 30 Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',...
    prob.getPercent(2,1,30,x),prob.getPercentile(1,0.95,x),t);
% OAR % > 30 Gy: 29.61, PTV D95: 79.03, Time: 33.43 (total 33.61)

% Load approximate dose with continuation (b)
load('ex1Results/ex1ContinueB.mat')
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose with continuation (b)')
fprintf('OAR %% > 30 Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',...
    prob.getPercent(2,1,30,x),prob.getPercentile(1,0.95,x),t);
% OAR % > 30 Gy: 19.84, PTV D95: 78.05, Time: 163.33 (total 163.51)

prob.plotDVHPaper(xMat)  
