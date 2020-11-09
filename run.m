% Example for setting up problem with FluenceMapOpt
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each organ included in the plan, create a structure containing the
% following fields:
% 
%   name: string used in data files
%   constraints: cell of organ constraints
%
% Each term has the following fields ('unif' doesn't need 'percent'):
% 
%   type: string 'unif', 'ldvc', or 'udvc'
%   dose: dose in Gy
%   percent: no more than p% receives more than (udvc) or less than (ldvc)
%       specified dose value
%   weight: weight in objective function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PTV - prostate
ptv.name = 'PTV_68';
ptv.terms = {struct('type','unif','dose',81,'weight',1)};

% OAR - rectum
oar.name = 'Rectum';
oar.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};

% Create problem instance
structs = {ptv,oar};
prob = FluenceMapOpt(structs);

% Calculate beamlet intensities
prob.calcBeams();

% Plot objective function
prob.plotObj();

% Plot dose-volume histogram
prob.plotDVH();

% Plot beamlet intensities
prob.plotBeams();

% Plot dose
prob.plotDose();
