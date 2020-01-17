% Example for setting up problem with FluenceMapOpt
clear all; close all; clc;

% Add path to data, CT images, and solver
addpath(genpath('PROSTATE'));
addpath(genpath('Prostate_Dicom'));
addpath(genpath('minConf'));

%% Problem setup

% Prostate
prostate.name = 'PTV_68';
pt1.type = 'unif'; pt1.dose = 81; pt1.weight = 1;
pt2.type = 'ldvc'; pt2.dose = 0.95*81; pt2.percent = 95; pt2.weight = 1;
pt3.type = 'udvc'; pt3.dose = 1.05*81; pt3.percent = 5; pt3.weight = 1;

% Lymph nodes
nodes.name = 'PTV_56';
nt1.type = 'unif'; nt1.dose = 60; nt1.weight = 1;
nt2.type = 'udvc'; nt2.dose = 65; nt2.percent = 0; nt2.weight = 100;

% Rectum
rectum.name = 'Rectum';
rt1.type = 'udvc'; rt1.dose = 50; rt1.percent = 50; rt1.weight = 1;
rt2.type = 'udvc'; rt2.dose = 30; rt2.percent = 30; rt2.weight = 1;
rt3.type = 'udvc'; rt3.dose = 10; rt3.percent = 10; rt3.weight = 1;
rt4.type = 'udvc'; rt4.dose = 60; rt4.percent = 20; rt4.weight = 1;
rt5.type = 'udvc'; rt5.dose = 20; rt5.percent = 60; rt5.weight = 1;
rt6.type = 'udvc'; rt6.dose = 75; rt6.percent = 0; rt6.weight = 1;

% Bladder
bladder.name = 'Bladder';
bt1.type = 'udvc'; bt1.dose = 30; bt1.percent = 30; bt1.weight = 1;

% Problem parameters
pars.angles = 0:52:358;
pars.lambda = 1e-8;
pars.maxIter = 500;
pars.overlap = false;
pars.tol = 1e-3;

%% Problem parameters
prostate.terms = {pt1,pt2,pt3};
nodes.terms = {nt1};
rectum.terms = {rt1};
bladder.terms = {bt1};
pars.structs = {prostate,rectum};

% Create problem instance
f = FluenceMapOpt(pars);

% Constraint generation method
f.constGen();
f.plotDVH();
