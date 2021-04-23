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

% OAR - right femoral head
femR.name = 'Rt_femoral_head';
femR.terms = {struct('type','udvc','dose',30,'percent',0,'weight',1)};

% Create problem instance
fprintf('\nExample 6\n');
structs = {prostate,nodes,rectum,bladder,femL,femR};
prob = FluenceMapOpt(structs,'tol',1e-2);

% Calculate approximate doses
fprintf('\nCalculating approximate doses\n');
labels = ['a' 'b' 'c'];
for ii = 1:length(labels)
    fprintf('\nExample 6%s\n\n',labels(ii));
    if ii == 1
        prob.calcBeams();
    elseif ii == 2
        prob.calcBeamsIter();
    else
        prob.calcBeamsSlack();
    end
    fprintf('\nIterations: %d, Time: %.2f\n',prob.nIter,prob.time);
    prob.saveResults(['ex6' labels(ii) 'Approx.mat']);
end

% Calculate polished doses
fprintf('\nCalculating polished doses\n');
for ii = 1:length(labels)
    fprintf('\nExample 6%s\n',labels(ii));
    if ii == 1
        load('ex6aApprox.mat');
        x = results.x;
        prob.calcBeamsPolish(x);
        t1 = prob.time;
    elseif ii == 2
        x0 = prob.x0;
        prob.calcBeamsPolish(x0);
        t1 = prob.time;
    else
        prob.calcBeamsConvex();
        t1 = prob.time;
        prob.calcBeamsPolish(prob.x);
        t1 = t1 + prob.time;
    end
    fprintf('\nObjective: %.4e, Time: %.2f\n',prob.getObj('unif'),t1);
    prob.saveResults(['ex6' labels(ii) 'Polish.mat']);
end

% Calculating approximate dose with continuation
fprintf('\nCalculating approximate dose with continuation\n')
prob = calcBeamsContinue(prob,structs,0.99,0.01,100,true,0,false);
fprintf('Iterations: %d, Time: %.2f\n',prob.nIter,prob.time);
prob.saveResults('ex6Continue.mat');
