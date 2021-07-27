% Figure 4: Organ contours for CORT dataset

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% Create problem instance

% prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% lymph nodes
nodes.name = 'PTV_56';
nodes.terms = {struct('type','unif','dose',60,'weight',1)};

% rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};

% bladder
bladder.name = 'Bladder';
bladder.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};

% left femoral head
lfem.name = 'Lt_femoral_head';
lfem.terms = {struct('type','udvc','dose',60,'percent',0,'weight',1)};

% right femoral head
rfem.name = 'Rt_femoral_head';
rfem.terms = {struct('type','udvc','dose',60,'percent',0,'weight',1)};

% Create problem instance
structs = {prostate,nodes,rectum,bladder,lfem,rfem};
prob = FluenceMapOpt(structs,'x0',zeros(986,1));

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
names = prob.getNames(prob.structs,prob.nStructs);
mask = prob.getMaskStruct(names,prob.overlap);
body50 = mask{end}(idx1,idx2,50);
imagesc(CT50), hold on

% Indexes
organs = [3 5 6 4 2 1 7];
colors = [2 7 7 3 4 1 5];

% Plot organ contours
for i = 1:length(mask)
   contour(mask{organs(i)}(idx1,idx2,50),1,...
       'Color',myLines(colors(7),:),'LineWidth',2); % colors(i)
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
