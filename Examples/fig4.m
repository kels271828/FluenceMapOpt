% Figure 4: Objective function contours of the nonconvex relaxation applied
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
f = FluenceMapOpt(pars);

%% Plot objective function contours for relaxed problem

for i = 2:-1:1
    
    if i == 1
       f.xInit = [0; 2e3];
       f.maxIter = 37;
       wLim = [-30 20];
       wStep = 10;
    else 
       f.maxIter = 26;
       wLim = [-5 20];
       wStep = 5;
    end
    
    plotContourW(f,wLim,wStep);
    
end
