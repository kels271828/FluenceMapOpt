% Figure 9: Dose-volume histograms for prostate tumor with uniform dose
% target of 81 Gy with a 30%, 30 Gy dose-volume constraint on the rectum
% and with a 10%, 10 Gy dose-volume constraint on the rectum.

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

labels = ['b' 'c'];
doses = [30 10];
for ii = 1:length(doses)
    fprintf('%d %%, %d Gy example\n\n',doses(ii),doses(ii));

    % PTV - prostate
    prostate.name = 'PTV_68';
    prostate.terms = {struct('type','unif','dose',81,'weight',1)};

    % OAR - rectum
    rectum.name = 'Rectum';
    rectum.terms = {struct('type','udvc','dose',doses(ii),...
        'percent',doses(ii),'weight',1)};

    % Create problem instance
    structs = {prostate,rectum};
    prob = FluenceMapOpt(structs);

    % Load approximate dose
    load(['ex1Results/ex1' labels(ii) '_x1.mat'])
    x0 = results.x0;
    x1 = results.x;
    disp('Approximate dose')
    fprintf('OAR %% > %d Gy: %.2f, PTV D95: %.2f\n\n',doses(ii),...
        prob.getPercent(2,1,x1),prob.getPercentile(prob.structs{1}.A*x1,0.95));
    % OAR % > 30 Gy: 33.37, PTV D95: 79.18
    % OAR % > 10 Gy: 23.85, PTV D95: 74.81

    % Load polished dose
    load(['ex1Results/ex1' labels(ii) '_x2.mat'])
    x2 = results.x;
    disp('Polished Dose')
    fprintf('OAR %% > 50 Gy: %.2f, PTV D95: %.2f\n\n',...
        prob.getPercent(2,1,x2),prob.getPercentile(prob.structs{1}.A*x2,0.95));
    % OAR % > 50 Gy: 29.98, PTV D95: 79.12
    % OAR % > 50 Gy: 8.19, PTV D95: 57.29

    % Load approximate dose with continuation
    load(['ex1Results/ex1' labels(ii) '_x3.mat'])
    x3 = results.x;
    disp('Approximate dose with continuation')
    fprintf('OAR %% > %d Gy: %.2f, PTV D95: %.2f\n\n',doses(ii),...
        prob.getPercent(2,1,x3),prob.getPercentile(prob.structs{1}.A*x3,0.95));
    % OAR % > 30 Gy: 30.95, PTV D95: 79.13
    % OAR % > 10 Gy: 11.35, PTV D95: 68.62

    % Plot dose-volume histograms
    prob.plotDVHPaper([x0 x1 x2 x3])
end
