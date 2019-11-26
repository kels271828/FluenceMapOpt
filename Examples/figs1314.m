% Figures 13 and 14: Description.

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% Calculate dose

% Prostate
prostate.name = 'PTV_68';
tt1.type = 'unif'; tt1.dose = 81; tt1.weight = 1; 
tt2.type = 'ldvc'; tt2.dose = 81; tt2.percent = 95; tt2.weight = 100;
prostate.terms = {tt1};

% Lymph nodes
nodes.name = 'PTV_56';
nt1.type = 'unif'; nt1.dose = 60; nt1.weight = 1;
nt2.type = 'udvc'; nt2.dose = 65; nt2.percent = 0; nt2.weight = 100;
nodes.terms = {nt1};

% Rectum
rectum.name = 'Rectum';
rt1.type = 'udvc'; rt1.dose = 50; rt1.percent = 50; rt1.weight = 1;
rectum.terms = {rt1};

% Bladder
bladder.name = 'Bladder';
bt1.type = 'udvc'; bt1.dose = 30; bt1.percent = 30; bt1.weight = 1;
bladder.terms = {bt1};

% Calculate beamlet intensities
pars.structs = {prostate,rectum,bladder,nodes};
pars.tol = 5e-5;
f = FluenceMapOpt(pars);
sum(f.structs{2}.A*f.x > 50)/f.structs{2}.nVoxels
sum(f.structs{3}.A*f.x > 30)/f.structs{3}.nVoxels
f.calcDose();

%% Plot DVH

f.plotDVHPaper();
sum(f.structs{2}.A*f.x > 50)/f.structs{2}.nVoxels
sum(f.structs{3}.A*f.x > 30)/f.structs{3}.nVoxels

%% Plot dose

f.plotDosePaper();

% Remove extra whitespace
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
width = outerpos(3) - ti(1) - ti(3);
height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom width height];

%% Plot beamlets

f.plotBeamletsPaper()
