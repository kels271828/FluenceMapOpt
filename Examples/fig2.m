% Figure 2: Objective function contours of the nonconvex relaxation applied
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

% Create new problem instance
pars.structs = {tumor,rectum};
pars.angles = [17,353];
pars.lambda = 5e-6;
pars.xInit = zeros(2,1);
f = FluenceMapOpt(pars);

%% Plot objective function contours for relaxed problem 

maxIterVals = [1 3 32];
    
% Relaxed problems (initial, middle, last iterate)
for i = 1:size(maxIterVals,2)

    % Calculate w values and plot objective function contours
    f.maxIter = maxIterVals(i);
    f.calcDose();

    % Plot region
    f.getd;
    x = f.A\f.d;
    plotContourX(pars.lambda,rt1.weight,f.structs{2}.terms{1}.w,x);
    
end
