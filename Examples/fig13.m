% Example 4: Comparisons with one PTV and one OAR

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

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};

% Create problem instance
structs = {prostate,rectum};
prob = FluenceMapOpt(structs,'tol',1e-2);

% Dose-volume objectives
fprintf('\nExample 4 with Dose-volume Objectives\n');
labels = ['a' 'b' 'c'];
xMat = [];
for ii = 1:length(labels)
    fprintf('\nExample 4%s\n',labels(ii));
    load(['ex4Results/ex4' labels(ii) 'Approx.mat'])
    x = results.x;
    t = results.time;
    n = results.nIter;
    xMat = [xMat x];
    fprintf('Rectum %% > 30 Gy: %.2f\n',prob.getPercent(2,1,30,x))
    fprintf('Prostate D95: %.2f, Prostate %% < 81 Gy: %.2f\n',...
        prob.getPercentile(1,0.95,x),prob.getPercent(1,2,81,x));
    fprintf('Iter: %d, Time: %.2f\n',n,t);
    
    % Example 4a
    % Rectum % > 30 Gy: 38.35
    % Prostate D95: 79.65, Prostate % < 81 Gy: 55.32
    % Iter: 5, Time: 1.65 (total 1.83)
    
    % Example 4b
    % Rectum % > 30 Gy: 56.55
    % Prostate D95: 78.74, Prostate % < 81 Gy: 44.33
    % Iter: 173, Time: 5.02 (total 5.20)
    
    % Example 4c
    % Rectum % > 30 Gy: 20.63
    % Prostate D95: 77.89, Prostate % < 81 Gy: 43.74
    % Iter: 5, Time: 198.70
end

% Plot dose-volume histograms
xMat = [xMat(:,2),xMat(:,1),xMat(:,3)];
prob.plotDVHPaper(xMat,false)

% Dose-volume constraints
fprintf('\nExample 4 with Dose-volume Constraints\n');
labels = ['a' 'b' 'c' 'd'];
xMat = [];
for ii = 1:length(labels)
    fprintf('\nExample 4%s\n',labels(ii));
    if ii == 4
        load('ex4Results/ex4Continue.mat')
    else
        load(['ex4Results/ex4' labels(ii) 'Polish.mat'])
    end
    x = results.x;
    prob.x = x;
    t = results.time;
    xMat = [xMat x];
    fprintf('Rectum %% > 30 Gy: %.2f\n',prob.getPercent(2,1,30,x))
    fprintf('Prostate D95: %.2f, Prostate %% < 81 Gy: %.2f\n',...
        prob.getPercentile(1,0.95,x),prob.getPercent(1,2,81,x));
    fprintf('Obj: %.2f, Time: %.2f\n',prob.getObj('unif'),t);
    
    % Example 4a
    % Rectum % > 30 Gy: 29.43
    % Prostate D95: 81.00, Prostate % < 81 Gy: 4.37
    % Obj: 6.95, Time: 63.74 (total 65.56)

    % Example 4b
    % Rectum % > 30 Gy: 28.40
    % Prostate D95: 81.00, Prostate % < 81 Gy: 3.26
    % Obj: 9.12, Time: 63.42 (total 63.60)

    % Example 4c
    % Rectum % > 30 Gy: 28.16
    % Prostate D95: 81.00, Prostate % < 81 Gy: 2.20
    % Obj: 8.94, Time: 57.82 (total 379.46)

    % Example 4d
    % Rectum % > 30 Gy: 29.43
    % Prostate D95: 81.32, Prostate % < 81 Gy: 2.97
    % Obj: 7.80, Time: 46.76 (total 46.94)
end

% Plot dose-volume histograms
xMat = [xMat(:,1),xMat(:,4),xMat(:,2),xMat(:,3)];
prob.plotDVHPaper(xMat,false)
