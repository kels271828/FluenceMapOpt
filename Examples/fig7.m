% Figure 7: OAR Dose-Volume Histograms

clear all; close all; clc;

% add data to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% OAR doesn't meet dose-volume constraint

% PTV
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% OAR
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};

% Calculate beamlets
structs = {prostate,rectum};
prob1 = FluenceMapOpt(structs);
prob1.plotDVHPaper();

% Voxels over 50 Gy
fprintf('Percent voxels > 50 Gy: %.2f\n',prob1.getPercent(2,1))

%% OAR does meet dose-volume constraint

% Adjust dose-volume constraint
rectum.terms = {struct('type','udvc','dose',45,'percent',45,'weight',1)};
structs = {prostate,rectum};
prob2 = FluenceMapOpt(structs);

% Calculate beamlets
prob2.calcBeams();
prob1.x0 = prob2.x;
prob1.x = prob2.x;
prob1.plotDVHPaper();

% Voxels over 50 Gy
fprintf('Percent voxels > 50 Gy: %.2f\n',prob1.getPercent(2,1))
