% Figure 12: Dose and beamlets for multiple PTVs and OARs

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% PTV - lymph nodes
nodes.name = 'PTV_56';
nodes.terms = {struct('type','unif','dose',60,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};

% OAR - bladder
bladder.name = 'Bladder';
bladder.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};

% Create problem instance
structs = {prostate,nodes,rectum,bladder};
prob = FluenceMapOpt(structs);

% Load approximate dose
load('ex3Results/ex3Approx.mat')
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
