% Figure 6: PTV Dose-Volume Histograms

clear all; close all; clc;

% add data to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% OAR doesn't meet dose-volume constraint

% PTV
tumor.name = 'PTV_68';
tt1.type = 'unif'; tt1.dose = 81; tt1.weight = 1; 
tumor.terms = {tt1};

% OAR
rectum.name = 'Rectum';
rt1.type = 'udvc'; rt1.dose = 50; rt1.percent = 50; rt1.weight = 1;
rectum.terms = {rt1};

% Calculate beamlets
pars.structs = {tumor,rectum};
f = FluenceMapOpt(pars);
f.plotDVHPaper();

% Voxels over 50 Gy
100*sum(f.structs{2}.A*f.x > 50)/f.structs{2}.nVoxels

%% OAR does meet dose-volume constraint

% Adjust dose-volume constraint
rt1.dose = 45; rt1.percent = 45;
rectum.terms = {rt1};
pars.structs = {tumor,rectum};
f = FluenceMapOpt(pars);

% Calculate beamlets
f.calcDose();
f.structs{2}.terms{1}.dose = 50;
f.structs{2}.terms{1}.percent = 50;
f.xInit = f.x;
f.plotDVHPaper();

% Voxels over 50 Gy
100*sum(f.structs{2}.A*f.x > 50)/f.structs{2}.nVoxels
