% Figure 4: Objective function contours of the nonconvex relaxation applied
% to the example in the introduction. 

% NOTE: Need to change SetAccess of FluenceMap Opt property structs.

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% Set up problem

% PTV - prostate voxel 1675228
prostate = struct('name','tumorEx');
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% OAR - rectum voxels 1674687 and 1675607
rectum = struct('name','rectumEx');
rectum.terms = {struct('type','udvc','dose',20,'percent',50,'weight',10)};

% Create problem instance
structs = {prostate,rectum};
prob = FluenceMapOpt(structs,'angles',[17 353],'lambda',5e-6,...
    'maxIter',32);

%% Plot objective function contours for relaxed problem

for i = 2:-1:1
    
    if i == 1
       prob.x0 = [0; 2e3];
       prob.maxIter = 37;
       wLim = [-30 20];
       wStep = 10;
    else 
       prob.maxIter = 26;
       wLim = [-5 20];
       wStep = 5;
    end
    
    plotContourW(prob,wLim,wStep);
    
end
