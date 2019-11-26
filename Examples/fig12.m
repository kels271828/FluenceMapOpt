% Figure 12: Prostate tumor with uniform dose target of 81 Gy with
% multiple dose-volume constraints.

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV
tumor.name = 'PTV_68';
tt1.type = 'unif'; tt1.dose = 81; tt1.weight = 1;
tumor.terms = {tt1};

% OAR
rectum.name = 'Rectum';
rt1.type = 'udvc'; rt1.dose = 20; rt1.percent = 60; rt1.weight = 1;
rt2.type = 'udvc'; rt2.dose = 50; rt2.percent = 50; rt2.weight = 1;
rt3.type = 'udvc'; rt3.dose = 60; rt3.percent = 20; rt3.weight = 1;
rt4.type = 'udvc'; rt4.dose = 75; rt4.percent = 0; rt4.weight = 1;
rectum.terms = {rt1,rt2,rt3,rt4};

% Initialize problem
pars.structs = {tumor,rectum};
pars.tol = 5e-5;
pars.maxIter = 1000;
f = FluenceMapOpt(pars);
sum(f.structs{2}.A*f.x > 20)/f.structs{2}.nVoxels
sum(f.structs{2}.A*f.x > 50)/f.structs{2}.nVoxels
sum(f.structs{2}.A*f.x > 60)/f.structs{2}.nVoxels
sum(f.structs{2}.A*f.x > 75)/f.structs{2}.nVoxels
f.calcDose();

% Plot dvh
f.plotDVHPaper();
sum(f.structs{2}.A*f.x > 20)/f.structs{2}.nVoxels
sum(f.structs{2}.A*f.x > 50)/f.structs{2}.nVoxels
sum(f.structs{2}.A*f.x > 60)/f.structs{2}.nVoxels
sum(f.structs{2}.A*f.x > 75)/f.structs{2}.nVoxels
