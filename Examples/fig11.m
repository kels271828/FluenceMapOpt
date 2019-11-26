% Figure 11: Prostate tumor with uniform dose target of 81 Gy with various
% dose-volume constraints on the rectum.

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% 30/30 DVC

% PTV
tumor.name = 'PTV_68';
tt1.type = 'unif'; tt1.dose = 81; tt1.weight = 1;
tumor.terms = {tt1};

% OAR
rectum.name = 'Rectum';
rt1.type = 'udvc'; rt1.dose = 30; rt1.percent = 30; rt1.weight = 1;
rectum.terms = {rt1};

% Set up problem
pars.tol = 5e-5;
pars.structs = {tumor,rectum};
f = FluenceMapOpt(pars);
sum(f.structs{2}.A*f.x > 30)/f.structs{2}.nVoxels

% Plot dose-volume histograms
f.calcDose();
f.plotDVHPaper();
sum(f.structs{2}.A*f.x > 30)/f.structs{2}.nVoxels

%% 10/10 DVC

% OAR
rt1.dose = 10; rt1.percent = 10;
rectum.terms = {rt1};

% Set up problem
pars.maxIter = 1000;
pars.structs = {tumor,rectum};
f = FluenceMapOpt(pars);
sum(f.structs{2}.A*f.x > 10)/f.structs{2}.nVoxels

% Plot dose-volume histograms
f.calcDose()
f.plotDVHPaper();
sum(f.structs{2}.A*f.x > 10)/f.structs{2}.nVoxels
