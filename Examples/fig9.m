% Figures 9-11: Dose-volume histograms for one PTV and one OAR with one 
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
fprintf('\nExample 1\n');
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
% OAR % > 30 Gy: 34.16, PTV D95: 79.17, Time: 6.38

% Load approximate dose with continuation 1
load('ex1Results/ex1Continue1.mat')
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose with continuation 1')
fprintf('OAR %% > 30 Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',...
    prob.getPercent(2,1,30,x),prob.getPercentile(1,0.95,x),t);
% OAR % > 30 Gy: 29.61, PTV D95: 79.03, Time: 0.84 (total 32.81)

% Load approximate dose with continuation 2
load('ex1Results/ex1Continue2.mat')
x = results.x;
xMat = [xMat x];
t = results.time;
disp('Approximate dose with continuation 2')
fprintf('OAR %% > 30 Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',...
    prob.getPercent(2,1,30,x),prob.getPercentile(1,0.95,x),t);
% OAR % > 30 Gy: 19.36, PTV D95: 78.02, Time: 3.37 (total 163.87)

prob.plotDVHPaper(xMat)  
