% Figure 5: Organ contours for CORT dataset

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% Create problem instance

% prostate
prostate.name = 'PTV_68';
tt1.type = 'unif'; tt1.dose = 81; tt1.weight = 1; 
prostate.terms = {tt1};

% lymph nodes
nodes.name = 'PTV_56';
nt1.type = 'unif'; nt1.dose = 60; nt1.weight = 1;
nodes.terms = {nt1};

% rectum
rectum.name = 'Rectum';
rt1.type = 'udvc'; rt1.dose = 50; rt1.percent = 50; rt1.weight = 1;
rectum.terms = {rt1};

% bladder
bladder.name = 'Bladder';
bt1.type = 'udvc'; bt1.dose = 30; bt1.percent = 30; bt1.weight = 1;
bladder.terms = {bt1};

% left femoral head
lfem.name = 'Lt_femoral_head';
lt1.type = 'udvc'; lt1.dose = 10; lt1.percent = 10; lt1.weight = 1;
lfem.terms = {lt1};

% right femoral head
rfem.name = 'Rt_femoral_head';
qt1.type = 'udvc'; qt1.dose = 10; qt1.percent = 10; qt1.weight = 1;
rfem.terms = {qt1};

% Create problem instance
pars.structs = {prostate,nodes,rectum,bladder,lfem,rfem};
f = FluenceMapOpt(pars);
f.x = zeros(size(f.x));

%% Plot structures

figure()
myLines = lines;
myLines(5,:) = [0 0 0];
idx1 = 40:126;
idx2 = 23:152;

% Get CT slice
ct = dicomread('CT.2.16.840.1.113662.2.12.0.3173.1271873797.276');
ct = double(imresize(ct,[184,184]));
ct50 = ct(idx1,idx2);
ctShift = ct50 - min(ct50(:));
ctShiftScale = ctShift/max(ctShift(:));
CT50 = repmat(ctShiftScale,[1 1 3]);

% Plot CT
body50 = f.mask{end}(idx1,idx2,50);
imagesc(CT50), hold on

% Indexes
organs = [3 5 6 4 2 1];
colors = [2 7 7 3 4 1];

% Plot organ contours
for i = 1:length(f.mask)-1
   contour(f.mask{organs(i)}(idx1,idx2,50),1,'Color',myLines(colors(i),:),'LineWidth',2); 
end

% Annotations
axis equal
axis off

% Remove extra whitespace
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
width = outerpos(3) - ti(1) - ti(3);
height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom width height];
