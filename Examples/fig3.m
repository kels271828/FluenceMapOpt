% Figure 3: Objective function contours of the nonconvex relaxation applied
% to the example in the introduction. 

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% Set up problem

% PTV - tumor voxel 1675228
tumor.name = 'tumorEx';
tt1.type = 'unif'; tt1.dose = 81; tt1.weight = 1;
tumor.terms = {tt1};

% OAR - rectum voxels 1674687 and 1675607
rectum.name = 'rectumEx';
rt1.type = 'udvc'; rt1.dose = 20; rt1.percent = 50; rt1.weight = 10;
rectum.terms = {rt1};

% Create problem instance
pars.structs = {tumor,rectum};
pars.angles = [17,353];
pars.lambda = 5e-6;
pars.maxIter = 32;
pars.xInit = zeros(2,1);
f = FluenceMapOpt(pars);

%% Plot objective function contours for relaxed problem

plotContourW(f,[-30 20],10);
