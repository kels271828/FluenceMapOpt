% Figure 8: Prostate tumor with uniform dose target of 81 Gy with a 50%,
% 50 Gy dose-volume constraint on the rectum.

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% Calculate dose

pars.tol = 5e-5;
f = FluenceMapOpt(pars);
f.calcDose();

%% Plot DVH

f.plotDVHPaper();
sum(f.structs{2}.A*f.x > 50)/f.structs{2}.nVoxels

%% Plot dose

f.plotDosePaper();

% Remove extra whitespace
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
width = outerpos(3) - ti(1) - ti(3);
height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom width height];

%% Plot beamlets

f.plotBeamletsPaper()

%% Plot objective values

f.plotObjPaper();
