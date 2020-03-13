% Figure 1: Example of a nonconvex feasible beamlet region due to
% dose-volume constraints. 

% To visualize the feasible region in 2D, we choose two beamlets (beam 
% angle 16 degrees, beamlet at index 86, and beam angle 352 degrees, 
% beamlet at index 85). The PTV target is a uniform dose of 81 Gy for tumor
% voxel 1675228, with a 20 Gy / 50% dose-volume constraint on rectum voxels
% 1674687 and 1675607.

% NOTE: Need to change SetAccess of FluenceMapOpt properties D and mask.

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV and OAR voxels
voxels = [1675228; 1674687; 1675607];

%% Find beamlets

% For tumor voxel 1675228 and rectum voxels 1674687, 1675607, figure out
% which two beamlets and angles are best for our examples
% vals = zeros(180,3);

% For each beam, find the beamlets that deposit the most radiation into the
% three voxels chosen for the example
% for i = 1:180
%     
%     angle = 2*(i-1);
%     load(['Gantry' int2str(angle) '_Couch0_D.mat']);
%     D = D(voxels,:);
%     [val,idx] = max(sum(D));
%     vals(i,:) = [angle, idx, full(val)];
% 
% end

% Beamlets chosen for example: 
%    16.0000   86.0000    0.0068
%   352.0000   85.0000    0.0064
%
% Beam angle 16 degrees, beamlet at index 86
% Beam angle 352 degrees, beamlet at index 85

%% Set up problem

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',80,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',20,'percent',50,'weight',1)};

% Creat problem instance
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);

%% Plot dose on slice 50

% Subsets of dose matrix
load('Gantry16_Couch0_D.mat');
temp = D(:,86);
A = D(voxels,86);
load('Gantry352_Couch0_D.mat');
prob.D = [temp D(:,85)];
A = [A D(voxels,85)];

% Add voxels to contour mask
myVoxels = zeros(184,184,90);
myVoxels(voxels) = 1;
prob.mask = [myVoxels prob.mask(:)'];

% Beamlet intensities at global minimum
lambda = 5e-6;
B1 = (20*lambda*A(3,1) + (81*A(3,2) - 20*A(1,2))*(A(1,1)*A(3,2) - ...
    A(1,2)*A(3,1)))/((A(1,1)*A(3,2) - A(1,2)*A(3,1))^2 +...
    lambda*(A(3,1)^2 + A(3,2)^2));
B2 = (20 - A(3,1)*B1)/A(3,2);
prob.x = [full(B1); B2];

% Plot dose
prob.plotDosePaper();

% Remove extra whitespace
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
width = outerpos(3) - ti(1) - ti(3);
height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom width height];

%% Plot nonconvex feasible region

plotContourX(lambda);
