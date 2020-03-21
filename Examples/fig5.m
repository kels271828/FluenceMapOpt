% Figure 5: PTV Dose-Volume Histograms

clear all; close all; clc;

% add data to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

%% Idealized uniform PTV dose of 81 Gy

myLines = lines;

% Plot dvh
plot([81 81],[0 100],':','Color',[0.4 0.4 0.4],'LineWidth',2); hold on
plot([0 81],[100 100],'Color',myLines(1,:),'LineWidth',2);
plot([81 100],[0 0],'Color',myLines(1,:),'LineWidth',2)
plot(81,100,'o','MarkerEdgeColor',myLines(1,:),'MarkerFaceColor',myLines(1,:),'MarkerSize',10)
plot(81,0,'o','MarkerEdgeColor',myLines(1,:),'MarkerFaceColor',[1 1 1],'MarkerSize',10)

% Annotations
ax = gca;
ax.XLim = [0 100];
ax.YLim = [0 100];
ax.XTick = 0:20:100;
ax.YTick = 0:20:100;
ax.XTickLabel = {};
ax.YTickLabel = {};
ax.LineWidth = 2;
axis square

%% Initialization for uniform PTV dose of 81 Gy

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% Calculate beamlets
structs = {prostate};
prob = FluenceMapOpt(structs);
prob.plotDVHPaper(prob.x);
