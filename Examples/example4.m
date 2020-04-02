% Example 4: Comparisons with multiple PTVs and OARs

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1)};

% PTV - lymph nodes
nodes.name = 'PTV_56';
nodes.terms = {struct('type','unif','dose',60,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',50,'percent',50,'weight',1)};

% OAR - bladder
bladder.name = 'Bladder';
bladder.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};

% Create problem instance
fprintf('\nExample 4\n');
structs = {prostate,nodes,rectum,bladder};
prob = FluenceMapOpt(structs);

% Calculate approximate doses
fprintf('\nCalculating approximate doses\n');
labels = ['b' 'c'];
for ii = 1:length(labels)
   fprintf('\nExample 4%s\n\n',labels(ii));
   if ii == 1                % (b) Iterative method
       prob.tol = 1e-2;
       prob.calcBeamsIter();
   else                      % (c) Slack method
       prob.tol = 1e-3;
       prob.calcBeamsSlack();
   end
   fprintf('\nIterations: %d, Time: %.2f\n',prob.nIter,prob.time);
   prob.saveResults(['ex4' labels(ii) 'Approx.mat']);
end

% Calculate polished doses
fprintf('\nCalculating polished doses\n');
labels = ['a' 'b' 'c'];
for ii = 1:length(labels)
   fprintf('\nExample 4%s\n',labels(ii));
   if ii == 1                    % (a) Our method
       load('ex3Results/ex3Approx.mat');
       x = results.x;
       prob.calcBeamsPolish(x);
   elseif ii == 2                % (b) Constraint generation method
       x0 = prob.x0;
       prob.calcBeamsPolish(x0);
   else                          % (c) Convex method
       prob.calcBeamsConvex();
       t1 = prob.time;
       prob.calcBeamsPolish(prob.x);
       prob.time = t1 + prob.time;
   end
   fprintf('\nTime: %.2f\n',prob.time);
   prob.saveResults(['ex4' labels(ii) 'Polish.mat']);
end
