% Example 5: Comparisons with multiple PTVs and OARs

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1),...
    struct('type','ldvc','dose',81,'percent',5,'weight',1)};

% PTV - lymph nodes
nodes.name = 'PTV_56';
nodes.terms = {struct('type','unif','dose',60,'weight',1),...
    struct('type','ldvc','dose',60,'percent',5,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};

% OAR - bladder
bladder.name = 'Bladder';
bladder.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};

% Create problem instance
structs = {prostate,nodes,rectum,bladder};
prob = FluenceMapOpt(structs,'tol',1e-2);

% Dose-volume objectives
fprintf('\nExample 5 with Dose-volume Objectives\n');
labels = ['a' 'b' 'c'];
xMat = [];
for ii = 1:length(labels)
    fprintf('\nExample 5%s\n',labels(ii));
    load(['ex5Results/ex5' labels(ii) 'Approx.mat'])
    x = results.x;
    t = results.time;
    n = results.nIter;
    xMat = [xMat x];
    fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 30 Gy: %.2f\n',...
        prob.getPercent(3,1,50,x),prob.getPercent(4,1,30,x));
    fprintf('Prostate D95: %.2f, Lymph Nodes D95: %.2f\n',...
        prob.getPercentile(1,0.95,x),prob.getPercentile(2,0.95,x));
    fprintf('Prostate %% < 81 Gy: %.2f, Lymph Nodes %% < 60 Gy: %.2f\n',...
        prob.getPercent(1,2,81,x),prob.getPercent(2,2,60,x));
    fprintf('Iter: %d, Time: %.2f\n',n,t);
    
    % Example 5a
    % Rectum % > 50 Gy: 64.58, Bladder % > 30 Gy: 49.47
    % Prostate D95: 76.95, Lymph Nodes D95: 57.51
    % Prostate % < 81 Gy: 57.86, Lymph Nodes % < 60 Gy: 52.95
    % Iter: 3, Time: 2.43 (total 2.80)

    % Example 5b
    % Rectum % > 50 Gy: 77.48, Bladder % > 30 Gy: 44.45
    % Prostate D95: 76.57, Lymph Nodes D95: 56.87
    % Prostate % < 81 Gy: 42.07, Lymph Nodes % < 60 Gy: 41.18
    % Iter: 195, Time: 12.61 (total 12.98)

    % Example 5c
    % Rectum % > 50 Gy: 25.68, Bladder % > 30 Gy: 27.55
    % Prostate D95: 74.18, Lymph Nodes D95: 57.34
    % Prostate % < 81 Gy: 45.73, Lymph Nodes % < 60 Gy: 42.25
    % Iter: 3, Time: 2014.94
end

% Plot dose-volume histograms
xMat = [xMat(:,2),xMat(:,1),xMat(:,3)];
prob.plotDVHPaper(xMat,false)

% Dose-volume constraints
fprintf('\nExample 5 with Dose-volume Constraints\n');
labels = ['a' 'c' 'd'];
xMat = [];
for ii = 1:length(labels)
    fprintf('\nExample 5%s\n',labels(ii));
    if ii == 3
        load('ex5Results/ex5Continue.mat')
    else
        load(['ex5Results/ex5' labels(ii) 'Polish.mat'])
    end
    x = results.x;
    prob.x = x;
    t = results.time;
    xMat = [xMat x];
    fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 30 Gy: %.2f\n',...
        prob.getPercent(3,1,50,x),prob.getPercent(4,1,30,x));
    fprintf('Prostate D95: %.2f, Lymp Nodes D95: %.2f\n',...
        prob.getPercentile(1,0.95,x),prob.getPercentile(2,0.95,x));
    fprintf('Prostate %% < 81 Gy: %.2f, Lymph Nodes %% < 60 Gy: %.2f\n',...
        prob.getPercent(1,2,81,x),prob.getPercent(2,2,60,x));
    fprintf('Obj: %.2f, Time: %.2f\n',prob.getObj('unif'),t); 
    
    % Example 5a
    % Rectum % > 50 Gy: 44.31, Bladder % > 30 Gy: 28.78
    % Prostate D95: 81.00, Lymp Nodes D95: 60.01
    % Prostate % < 81 Gy: 4.46, Lymph Nodes % < 60 Gy: 3.68
    % Obj: 28.64, Time: 307.03 (total 309.83)

    % Example 5c
    % Rectum % > 50 Gy: 33.41, Bladder % > 30 Gy: 25.02
    % Prostate D95: 81.00, Lymp Nodes D95: 60.01
    % Prostate % < 81 Gy: 3.07, Lymph Nodes % < 60 Gy: 2.69
    % Obj: 36.02, Time: 193.61 (total 1118.13)

    % Example 5d
    % Rectum % > 50 Gy: 48.39, Bladder % > 30 Gy: 29.84
    % Prostate D95: 81.17, Lymp Nodes D95: 60.32
    % Prostate % < 81 Gy: 4.71, Lymph Nodes % < 60 Gy: 3.90
    % Obj: 27.71, Time: 137.99 (total 138.36)
end

% Plot dose-volume histograms
xMat = [xMat(:,1),xMat(:,3),xMat(:,2)];
prob.plotDVHPaper(xMat,false)