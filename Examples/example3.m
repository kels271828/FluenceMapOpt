% Example 3: Multiple PTVs and OARs

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
structs = {prostate,nodes,rectum,bladder};
prob = FluenceMapOpt(structs);

% Loop over different methods
labels = ['a' 'b' 'c' 'd'];
for ii = 1:4
    fprintf('\nExample 3%s\n',labels(ii));
    
    % Calculate approximate dose
    fprintf('\nCalculating approximate dose\n\n');
    if ii == 1
        prob.calcBeams();       % a) Our method    
    elseif ii == 2
        prob.calcBeamsConvex(); % b) Convex method
    elseif ii == 3
        prob.calcBeamsIter();   % c) Iterative method
    else
        prob.calcBeamsSlack();  % d) Slack method
    end
    
    % Approximate dose results
    if ii == 2
        fprintf('\nTime: %.2f\n',prob.time);
    else
        fprintf('\nIterations: %d, Time: %.2f\n',prob.nIter,prob.time);
    end
    filename1 = ['e3' labels(ii) 'Approx.mat'];
    prob.saveResults(filename1);
    
    % Calculate polished dose
    fprintf('\nCalculating polished dose\n');
    prob.calcBeamsPolish(prob.x);
    fprintf('\nTime: %.2f\n',prob.time);
    filename2 = ['e3' labels(ii) 'Polish.mat'];
    prob.saveResults(filename2);
end

% f) Polish initialization
fprintf('\nExample 3e\n');
fprintf('\nCalculating polished dose\n');
prob.calcBeamsPolish(prob.x0);
prob.saveResults('e3ePolish.mat');
