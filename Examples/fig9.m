% Figure 9: Dose-volume histograms for one PTV and one OAR with one 
% dose-volume constraint

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1),...
    struct('type','udvc','dose',30,'percent',30,'weight',1),...
    struct('type','udvc','dose',10,'percent',10,'weight',1)};

% Create problem instance
structs = {prostate,rectum};
prob = FluenceMapOpt(structs);
xMat = prob.x0;
disp('Initialization')
fprintf('OAR %% > 50 Gy: %.2f, %% > 30 Gy: %.2f, %% > 10 Gy: %.2f, PTV D95: %.2f\n\n',...
    prob.getPercent(2,1,xMat),prob.getPercent(2,2,xMat),prob.getPercent(2,3,xMat),...
    prob.getPercentile(prob.structs{1}.A*xMat,0.95));
% OAR % > 50 Gy: 56.80, % > 30 Gy: 64.14, % > 10 Gy: 73.97, PTV D95: 79.65

labels = ['a' 'b' 'c'];
doses = [50 30 10];
for ii = 1:length(labels)
    % Load approximate dose
    filename = ['ex1Results/ex1' labels(ii) 'Approx.mat'];
    load(filename)
    x = results.x;
    xMat = [xMat x];
    t = results.time;
    fprintf('\nExample 1%s\n',labels(ii));
    fprintf('OAR %% > %d Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',...
        doses(ii),prob.getPercent(2,ii,x),...
        prob.getPercentile(prob.structs{1}.A*x,0.95),t);
    % OAR % > 50 Gy: 51.52, PTV D95: 79.66, Time: 23.69
    % OAR % > 30 Gy: 33.37, PTV D95: 79.18, Time: 47.28
    % OAR % > 10 Gy: 23.85, PTV D95: 74.81, Time: 33.48
end

% Plot dose-volume histograms
prob.plotDVHPaper(xMat)
