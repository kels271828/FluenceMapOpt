% Example 6: Comparisons with multiple PTVs and OARs

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
bladder.terms = {struct('type','udvc','dose',40,'percent',40,'weight',1)};

% OAR - left femoral head
femL.name = 'Lt_femoral_head';
femL.terms = {struct('type','udvc','dose',30,'percent',0,'weight',1)};

femR.name = 'Rt_femoral_head';
femR.terms = {struct('type','udvc','dose',30,'percent',0,'weight',1)};

% Create problem instance
structs = {prostate,nodes,rectum,bladder,femL,femR};
prob = FluenceMapOpt(structs,'tol',1e-2);
% xInit = 0.298422 seconds

% Dose-volume objectives
fprintf('\nExample 6 with Dose-volume Objectives\n');
labels = ['a' 'b' 'c'];
xMat = [];
for ii = 1:length(labels)
    fprintf('\nExample 6%s\n',labels(ii));
    load(['ex6Results/ex6' labels(ii) 'Approx.mat'])
    x = results.x;
    t = results.time;
    n = results.nIter;
    xMat = [xMat x];
    fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 40 Gy: %.2f\n',...
        prob.getPercent(3,1,50,x),prob.getPercent(4,1,40,x));
    fprintf('LFemur %% > 30 Gy: %.2f, RFemur %% > 30 Gy: %.2f\n',...
        prob.getPercent(5,1,30,x),prob.getPercent(6,1,30,x));
    fprintf('Prostate D95: %.2f, Lymph Nodes D95: %.2f\n',...
        prob.getPercentile(1,0.95,x),prob.getPercentile(2,0.95,x));
    fprintf('Prostate %% < 81 Gy: %.2f, Lymph Nodes %% < 60 Gy: %.2f\n',...
        prob.getPercent(1,2,81,x),prob.getPercent(2,2,60,x));
    fprintf('Iter: %d, Time: %.2f\n',n,t);
    
    % Example 6a
    % Rectum % > 50 Gy: 65.61, Bladder % > 40 Gy: 52.89
    % LFemur % > 30 Gy: 2.37, RFemur % > 30 Gy: 2.38
    % Prostate D95: 76.83, Lymph Nodes D95: 58.25
    % Prostate % < 81 Gy: 58.18, Lymph Nodes % < 60 Gy: 53.77
    % Iter: 2, Time: 1.86 (total 2.16)
    % 
    % Example 6b
    % Rectum % > 50 Gy: 73.89, Bladder % > 40 Gy: 53.02
    % LFemur % > 30 Gy: 3.88, RFemur % > 30 Gy: 3.63
    % Prostate D95: 76.86, Lymph Nodes D95: 58.13
    % Prostate % < 81 Gy: 43.68, Lymph Nodes % < 60 Gy: 46.09
    % Iter: 97, Time: 7.08 (total 7.38)
    % 
    % Example 6c
    % Rectum % > 50 Gy: 21.06, Bladder % > 40 Gy: 19.37
    % LFemur % > 30 Gy: 1.65, RFemur % > 30 Gy: 1.79
    % Prostate D95: 74.17, Lymph Nodes D95: 57.78
    % Prostate % < 81 Gy: 45.21, Lymph Nodes % < 60 Gy: 42.80
    % Iter: 2, Time: 29246.23
end

% Plot dose-volume histograms
xMat = [xMat(:,2),xMat(:,1),xMat(:,3)];
prob.plotDVHPaper(xMat,false)

% Dose-volume constraints
fprintf('\nExample 6 with Dose-volume Constraints\n');
labels = ['a' 'b' 'd'];
xMat = [];
for ii = 1:length(labels)
    fprintf('\nExample 6%s\n',labels(ii));
    if ii == 3
        load('ex6Results/ex6Continue.mat')
    else
        load(['ex6Results/ex6' labels(ii) 'Polish.mat'])
    end
    x = results.x;
    prob.x = x;
    t = results.time;
    xMat = [xMat x];
    fprintf('Rectum %% > 50 Gy: %.2f, Bladder %% > 40 Gy: %.2f\n',...
        prob.getPercent(3,1,50,x),prob.getPercent(4,1,40,x));
    fprintf('LFemur %% > 30 Gy: %.2f, RFemur %% > 30 Gy: %.2f\n',...
        prob.getPercent(5,1,30,x),prob.getPercent(6,1,30,x));
    fprintf('Prostate D95: %.2f, Lymp Nodes D95: %.2f\n',...
        prob.getPercentile(1,0.95,x),prob.getPercentile(2,0.95,x));
    fprintf('Prostate %% < 81 Gy: %.2f, Lymph Nodes %% < 60 Gy: %.2f\n',...
        prob.getPercent(1,2,81,x),prob.getPercent(2,2,60,x));
    fprintf('Obj: %.2f, Time: %.2f\n',prob.getObj('unif'),t); 
    
    % Example 6a
    % Rectum % > 50 Gy: 48.27, Bladder % > 40 Gy: 37.58
    % LFemur % > 30 Gy: 0.00, RFemur % > 30 Gy: 0.00
    % Prostate D95: 81.00, Lymp Nodes D95: 60.01
    % Prostate % < 81 Gy: 4.89, Lymph Nodes % < 60 Gy: 4.10
    % Obj: 15.25, Time: 414.99 (total 417.34)
    % 
    % Example 6b
    % Rectum % > 50 Gy: 37.49, Bladder % > 40 Gy: 32.75
    % LFemur % > 30 Gy: 0.00, RFemur % > 30 Gy: 0.00
    % Prostate D95: 81.01, Lymp Nodes D95: 60.26
    % Prostate % < 81 Gy: 4.17, Lymph Nodes % < 60 Gy: 1.19
    % Obj: 52.90, Time: 473.96 (total 474.25)
    % 
    % Example 6d
    % Rectum % > 50 Gy: 46.68, Bladder % > 40 Gy: 38.12
    % LFemur % > 30 Gy: 0.00, RFemur % > 30 Gy: 0.00
    % Prostate D95: 81.05, Lymp Nodes D95: 60.18
    % Prostate % < 81 Gy: 4.93, Lymph Nodes % < 60 Gy: 3.97
    % Obj: 23.67, Time: 215.08 (total 215.38)
end

% Plot dose-volume histograms
xMat = [xMat(:,1),xMat(:,3),xMat(:,2)];
prob.plotDVHPaper(xMat,false)