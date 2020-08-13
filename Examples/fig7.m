% Figure 7: Dose and beamlets for prostate tumor with uniform dose target
% of 81 Gy with a 30%, 30 Gy dose-volume constraint on the rectum.

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
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);

% Load approximate dose
load('ex1Results/ex1Approx.mat')
prob.x = results.x;

% Plot dose
prob.plotDosePaper()

% Remove extra whitespace
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
width = outerpos(3) - ti(1) - ti(3);
height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom width height];

% Plot beamlets
prob.plotBeamsPaper()
