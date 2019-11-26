% Example for setting up problem with FluenceMapOpt
clear all; close all; clc;

% Add path to data, CT images, and solver
addpath(genpath('PROSTATE'));
addpath(genpath('Prostate_Dicom'));
addpath(genpath('minConf'));

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
%   percent: at least (ldvc) or at most (udvc) p% receives at least d Gy
%   weight: weight in objective function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To run the default problem, there is no need to include the parameter
% structure 'pars'; just use the command f = FluenceMapOpt(). Below we set
% up the default problem by explicitly defining the input parameters. 

% Planning target volume (prostate)
ptv.name = 'PTV_68';
pt1.type = 'unif'; pt1.dose = 81; pt1.weight = 1;
ptv.terms = {pt1};

% Organ-at-risk (rectum)
oar.name = 'Rectum';
rt1.type = 'udvc'; rt1.dose = 50; rt1.percent = 50; rt1.weight = 1;
oar.terms = {rt1};

% Problem parameters
pars.structs = {ptv,oar};
pars.angles = 0:52:358;
pars.lambda = 1e-8;
pars.maxIter = 500;
pars.overlap = false;
pars.tol = 1e-3;

% Create problem instance
f = FluenceMapOpt(pars);

% Calculate beamlet intensities
f.calcDose();

% Plot objective function
f.plotObj();

% Plot dose-volume histogram
f.plotDVH();

% Plot beamlet intensities
f.plotBeamlets();

% Plot dose
f.plotDose();
