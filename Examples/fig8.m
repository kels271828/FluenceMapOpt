% Figure 8: Dose-volume histograms for prostate tumor with uniform dose
% target of 81 Gy with a 50%, 50 Gy dose-volume constraint on the rectum.

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% OAR
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};

% Create problem instance
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);

% Load approximate dose
load('ex1Results/ex1a_x1.mat')
x0 = results.x0;
x1 = results.x;
disp('Approximate dose')
fprintf('OAR %% > 50 Gy: %.2f, PTV D95: %.2f\n\n',...
    prob.getPercent(2,1,x1),prob.getPercentile(prob.structs{1}.A*x1,0.95));
% OAR % > 50 Gy: 51.52, PTV D95: 79.66

% Load polished dose
load('ex1Results/ex1a_x2.mat')
x2 = results.x;
disp('Polished Dose')
fprintf('OAR %% > 50 Gy: %.2f, PTV D95: %.2f\n\n',...
    prob.getPercent(2,1,x2),prob.getPercentile(prob.structs{1}.A*x2,0.95));
% OAR % > 50 Gy: 50.00, PTV D95: 79.65

% Load approximate dose with continuation
load('ex1Results/ex1a_x3.mat')
x3 = results.x;
disp('Approximate dose with continuation')
fprintf('OAR %% > 50 Gy: %.2f, PTV D95: %.2f\n',...
    prob.getPercent(2,1,x3),prob.getPercentile(prob.structs{1}.A*x3,0.95));
% OAR % > 50 Gy: 50.49, PTV D95: 79.66

% Plot dose-volume histograms
prob.plotDVHPaper([x0 x1 x3 x2])
