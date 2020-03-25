% Figure 8: Dose-volume histogram for one PTV and one OAR with one 
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

labels = ['b' 'c'];
doses = [30 10];
for ii = 1:length(labels)
    fprintf('Example 1%s\n\n',labels(ii));
    
    % OAR - rectum
    rectum.name = 'Rectum';
    rectum.terms = {struct('type','udvc','dose',doses(ii),'percent',...
        doses(ii),'weight',1)};

    % Create problem instance
    structs = {prostate,rectum};
    prob = FluenceMapOpt(structs);
    xMat = prob.x0;
    disp('Initialization')
    fprintf('OAR %% > %d Gy: %.2f, PTV D95: %.2f\n\n',doses(ii),...
        prob.getPercent(2,1,xMat),prob.getPercentile(prob.structs{1}.A*xMat,0.95));
    % OAR % > 30 Gy: 64.14, PTV D95: 79.65, Time: 0.1792
    % OAR % > 10 Gy: 73.97, PTV D95: 79.65, Time: 0.1792
        
    % Load approximate dose
    load(['ex1Results/ex1' labels(ii) 'Approx.mat'])
    x = results.x;
    xMat = [xMat x];
    t = results.time;
    disp('Approximate dose')
    fprintf('OAR %% > %d Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',doses(ii),...
    prob.getPercent(2,1,x),prob.getPercentile(prob.structs{1}.A*x,0.95),t);
    % OAR % > 30 Gy: 34.16, PTV D95: 79.17, Time: 6.74
    % OAR % > 10 Gy: 24.76, PTV D95: 74.77, Time: 4.38
    
    % Load approximate dose with continuation
    load(['ex1Results/ex1' labels(ii) 'Continue.mat'])
    x = results.x;
    xMat = [xMat x];
    t = results.time;
    disp('Approximate dose with continuation')
    fprintf('OAR %% > %d Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',doses(ii),...
        prob.getPercent(2,1,x),prob.getPercentile(prob.structs{1}.A*x,0.95),t);
    % OAR % > 30 Gy: 31.07, PTV D95: 79.11, Time: 428.06
    % OAR % > 10 Gy: 11.29, PTV D95: 68.18, Time: 1263.68
    
    % Load polished dose
    load(['ex1Results/ex1' labels(ii) 'Polish.mat'])
    x = results.x;
    xMat = [xMat x];
    t = results.time;
    disp('Polished dose')
    fprintf('OAR %% > %d Gy: %.2f, PTV D95: %.2f, Time: %.2f\n\n',doses(ii),...
        prob.getPercent(2,1,x),prob.getPercentile(prob.structs{1}.A*x,0.95),t);
    % OAR % > 30 Gy: 29.98, PTV D95: 79.15, Time: 7.37
    % OAR % > 10 Gy: 8.19, PTV D95: 57.29, Time: 9.20

    % Plot dose-volume histograms
    prob.plotDVHPaper(xMat)  
end
