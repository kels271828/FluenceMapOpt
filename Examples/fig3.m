% Figure 3: Objective function contours of the nonconvex relaxation applied
% to the example in the introduction. 

% NOTE: Need to change SetAccess of FluenceMapOpt property structs.

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
prob = FluenceMapOpt(structs,'angles',[17 353],'lambda',5e-6,...
    'maxIter',32,'x0',zeros(2,1));

%% Plot objective function contours for relaxed problem

plotContourW(prob,[-30 20],10);
