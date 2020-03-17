% Figure 8: Dose and beamlets for prostate tumor with uniform dose target
% of 81 Gy with a 50%, 50 Gy dose-volume constraint on the rectum.

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

% Load approximate dose
load('ex1Results/ex1a_x1.mat')
x0 = results.x0;
x1 = results.x;
prob.x = x1;

% Plot dose
prob.plotDosePaper()

% Plot beamlets
prob.plotBeamsPaper()
