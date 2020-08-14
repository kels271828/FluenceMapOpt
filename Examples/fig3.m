% Figure 3: Objective function contours of the nonconvex relaxation applied
% to the example in the introduction. 

% NOTE: Need to add the following files to PROSTATE folder:
% * Gantry17_Couch0_D.mat
% * Gantry353_Couch0_D.mat
% * rectumEx_VOILIST.mat
% * tumorEx_VOILIST.mat

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% Set up problem

% PTV - prostate voxel 1675228
prostate.name = 'tumorEx';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% OAR - rectum voxels 1674687 and 1675607
rectum.name = 'rectumEx';
rectum.terms = {struct('type','udvc','dose',20,'percent',50,'weight',10)};

% Create problem instance
structs = {prostate,rectum};
prob = FluenceMapOpt(structs,'angles',[17 353],'lambda',5e-6);

%% Plot objective function contours for relaxed problem

for ii = 1:3
    if ii == 1
        prob.maxIter = 26;
        wLim = [-5 20];
        wStep = 5;
    elseif ii == 2
        prob.x0 = zeros(2,1);
        prob.maxIter = 32;
        wLim = [-30 20];
        wStep = 10;
    else
        prob.x0 = [0; 2e3];
        prob.maxIter = 37;
        wLim = [-30 20];
        wStep = 10;
    end
    plotContourW(prob,structs,wLim,wStep);
end
