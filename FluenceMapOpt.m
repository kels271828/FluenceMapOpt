classdef FluenceMapOpt < handle
% FLUENCEMAPOPT Fluence map optimization with dose-volume constraints.
%
%   Problem statement:
% 
%   min_(x,w) 
%       sum(i in I) weight_i/(2*nVoxels_i)*||A_i*x - d_i||_2^2
%       + sum(j in J) weight_j/(2*nVoxels_j)*||w_j - (A_j*x - d_j)||_2^2
%       + lambda/2*||x||_2^2
%   s.t. 
%       x >= 0
%       ||max(0,w_j)||_0 <= nVoxels_j*percent_j/100 for all j in J
%   where
%       I = set of uniform dose targets
%       J = set of dose-volume constraints
%
%   For each body structure included in the treatment plan, create a
%   structure with the following fields:
%
%       name: string used in data files
%       terms: cell containing organ constraint terms
%
%   Each term should have the following fields:
%
%       type: string 'unif', 'ldvc', or 'udvc'
%       weight: weight coefficient of the term in the objective function
%       dose: dose in Gy
%       percent: include if type is 'ldvc' or 'udvc'
%           
%           'ldvc': No more than p% receives less than d Gy
%           'udvc': No more than p% receives more than d Gy
%
%   Written to work with the CORT prostate tumor dataset, but could be
%   modified to work with other datasets. 

    properties
        % Input parameters
        structs            % Body structures
        xInit              % Initial beamlet intensities
        angles = 0:52:358; % Gantry angles
        lambda = 1e-8;     % L2 regularization coefficient
        maxIter = 500;     % Maximum number of iterations
        tol = 1e-3;        % Stopping tolerance
        overlap = false;   % Allow overlaps in structures
        
        % Internal variables
        D        % Full beamlet-to-voxel matrix
        names    % Body structure names
        mask     % Body structure contours
        A        % Stacked beamlet-to-voxel matrix
        d        % Stacked dose vector
        lb       % Lower bound for beamlets
        ub       % Upper bound for beamlets
        
        % might want to add back nStructs,nTerms,nAngles,nBeamlets for
        % readibility
        
        % Solution variables
        x        % Final beamlet intensities
%         obj      % Objective function values
%         err      % Error values between w vectors
%         nIter    % Number of iterations used
    end
    
    methods
        function prob = FluenceMapOpt(structs,varargin)
        % FLUENCEMAPOPT Initialize problem.
        %
        %   prob = FluenceMapOpt(structs)
        %       Initialize problem with default parameters.
        %
        %   prob = FluenceMapOpt(structs,xInit,angles,lambda,maxIter,tol,overlap)
        %       Initialize problem with optional arguments.
        
            % Set input variables
            if nargin == 0
                error('Not enough input arguments.')
            end
            if ~iscell(structs)
                error('Invalid input for `structs`.')
            end
            prob.setInputVars(varargin,nargin-1);
            
            % Comput internal variables
            prob.D = getD(prob.angles);
            prob.structs = getStructVars(structs,prob.overlap,prob.D);
            prob.names = getNames(prob.structs);
            prob.mask = getMaskStruct(prob.names,prob.overlap);
            [prob.A,prob.lb,prob.ub] = getA(prob.structs,prob.lambda);  
            
            % Compute initial beamlets
            if nargin > 1
                prob.xInit = varargin{1};
            else
                prob.xInit = getInitX(prob.structs,prob.lambda);
            end
        end
        
        function setInputVars(prob,args,nArgs)
        % SETINPUTVARS Set input variables.
        
            varNames = {'angles','lambda','maxIter','tol','overlap'};
            for ii = 1:length(varNames)
               if nArgs > ii
                   if ~isempty(args{ii+1})
                       prob.(varNames{ii}) = args{ii+1};
                   end
               end
            end
        end
        
        % in progress...
        function calcBeamlets(prob,varargin)
        % CALCBEAMLETS Calculate beamlet intensities.
        %
        %   calcBeamlets()
        %       Calculate beamlets and print iteration output.
        %
        %   calcBeamlets(print)
        %       Calculate beamlets with or without iteration output.
            
            if nargin > 1
                print = varargin{1};
            else
                print = true;
            end
            
            % Fluence map optimization
            prob.initProb(print);
           
%             for t = 1:prob.maxIter
% 
%                 % Update x
%                 prob.projX();
%                 
%                 % Update w
%                 errorSum = 0;
%                 for ii = 1:length(f.structs)
%                     for jj = 1:length(prob.structs{ii}.terms)
%                         if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
%                             
%                             % Update wPrev
%                             wPrev = prob.structs{ii}.terms{jj}.w;
%                             
%                             % Update w
%                             Axmd = prob.structs{ii}.A*f.x - f.structs{ii}.terms{jj}.d;
%                             coeff = f.structs{ii}.terms{jj}.step*f.structs{ii}.terms{jj}.weight/f.structs{ii}.nVoxels;
%                             temp = f.structs{ii}.terms{jj}.w + coeff*(Axmd - f.structs{ii}.terms{jj}.w);
%                             s = strcmp(f.structs{ii}.terms{jj}.type,'ldvc');
%                             f.structs{ii}.terms{jj}.w = (-1)^s*f.projW((-1)^s*temp,f.structs{ii}.terms{jj}.k);
%                             
%                             % Error sum
%                             errorSum = errorSum + norm(f.structs{ii}.terms{jj}.w - wPrev)/f.structs{ii}.terms{jj}.step;  
%                         end
%                     end
%                 end
%                 prob.err(t) = errorSum;
%                 prob.nIter = t;
%                 prob.calcObj(t,print);
% 
%                 % Stopping criteria
%                 if errorSum <= prob.tol
%                     break
%                 end
        end
        
        % in progress...
        function initProb(prob,print)
        % INITPROB Initialize x, w, and objective values.
        
            prob.x = prob.xInit;
            prob.initW();
            
%             % Initialize objective function values
%             f.obj = zeros(1,f.maxIter+1);
%             f.err = zeros(1,f.maxIter);
%             for i = 1:f.nStructs
%                 for j = 1:f.structs{i}.nTerms
%                     f.structs{i}.terms{j}.obj = zeros(1,f.maxIter+1);
%                     if ~strcmp(f.structs{i}.terms{j}.type,'unif')
%                         f.structs{i}.terms{j}.vdiff = zeros(1,f.maxIter+1);
%                         f.structs{i}.terms{j}.wdiff = zeros(1,f.maxIter+1);
%                     end
%                 end
%             end
%             
%             % Calculate and print initial objective value
%             f.calcObj(0,print)
        end
        
        % need to test...
        function initW(prob)
        % INITW Initialize w vectors for dose-volume constraint terms.
        
            for ii = 1:length(prob.structs)
               for jj = 1:length(prob.structs{ii}.terms)
                   if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                       initDose = prob.structs{ii}.A*f.x;
                       diff = initDose - prob.structs{ii}.terms{jj}.d;
                       s = strcmp(prob.structs{ii}.terms{jj}.type,'ldvc');
                       w = projW((-1)^s*diff,prob.structs{ii}.terms{jj}.k);
                       prob.structs{ii}.terms{jj}.w = (-1)^s*w;
                   end
               end
            end
        end
        
        % in progress...
        function calcObj(prob,iter,print)
        % CALCOBJ Calculate and print objective function value.
        
            for ii = 1:length(prob.structs)
                for jj = 1:length(prob.structs{ii}.terms)
                    dose = prob.structs{ii}.A*prob.x;
                    diff = dose - prob.structs{ii}.terms{jj}.d;
                    if strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        % do something
                    else
                        % do something else
                    end
                end
            end
            
            for i = 1:f.nStructs
                for j = 1:f.structs{i}.nTerms
                    Axmd = f.structs{i}.A*f.x - f.structs{i}.terms{j}.d;
                    if strcmp(f.structs{i}.terms{j}.type,'unif')
                        f.structs{i}.terms{j}.obj(iter+1) = f.structs{i}.terms{j}.weight*norm(Axmd)^2/(2*f.structs{i}.nVoxels);
                    else
                        s = strcmp(f.structs{i}.terms{j}.type,'ldvc');
                        f.structs{i}.terms{j}.vdiff(iter+1) = 100*sum((-1)^s*Axmd > 0)/f.structs{i}.nVoxels;
                        f.structs{i}.terms{j}.wdiff(iter+1) = 100*sum((-1)^s*f.structs{i}.terms{j}.w > 0)/f.structs{i}.nVoxels;
                        f.structs{i}.terms{j}.obj(iter+1) = f.structs{i}.terms{j}.weight*norm(Axmd - f.structs{i}.terms{j}.w)^2/(2*f.structs{i}.nVoxels);
                    end
                    f.obj(iter+1) = f.obj(iter+1) + f.structs{i}.terms{j}.obj(iter+1);
                end
            end
            f.obj(iter+1) = f.obj(iter+1) + 0.5*f.lambda*norm(f.x)^2;
            % f.getd('full'); f.obj(iter+1) = 0.5*norm(f.A*f.x - f.d)^2;
            if print
                fprintf('iter: %d, obj: %7.4e\n',iter,f.obj(iter+1));
            end
        end
        
        function notYetChecked(prob)
%         % Constraint generation method.
%         function constGen(f)
%             
%             
%             [A_dvc,d_dvc] = f.calcConstMatBig(f.xInit);
%             f.x = lsqlin(f.A_unif,f.d_unif,A_dvc,d_dvc,[],[],f.lb,f.ub);
%         end
%         
%         % Calculate stacked constraint matrix for dose-volume constraints.
%         function [A_dvc,d_dvc] = calcConstMatBig(f,x)
%            
%             A_dvc = [];
%             d_dvc = [];
%             for i = 1:f.nStructs
%                 for j = 1:f.structs{i}.nTerms
%                     if ~strcmp(f.structs{i}.terms{j}.type,'unif')
%                         [A_temp,d_temp] = calcConstMatSmall(f,x,i,j);
%                         A_dvc = [A_dvc; A_temp];
%                         d_dvc = [d_dvc; d_temp];
%                     end
%                 end
%             end
%         end
%         
%         % Calculate constraint matrix for dose-volume constraint.
%         function [A_dvc,d_dvc] = calcConstMatSmall(f,x,struct,term)
%             
%             [~,idx] = sort(f.structs{struct}.A*x);
%             if strcmp(f.structs{struct}.terms{term}.type,'udvc')
%                 idx_lower = f.structs{struct}.nVoxels - f.structs{struct}.terms{term}.k;
%                 A_dvc = f.structs{struct}.A(idx(1:idx_lower),:);
%                 d_dvc = f.structs{struct}.terms{term}.d(idx(1:idx_lower));
%             elseif strcmp(f.structs{struct}.terms{term}.type,'ldvc')
%                 idx_upper = f.structs{struct}.terms{term}.k + 1;
%                 A_dvc = -f.structs{struct}.A(idx(idx_upper:end),:);
%                 d_dvc = -f.structs{struct}.terms{term}.d(idx(idx_upper:end));
%             else
%                 A_dvc = NaN;
%                 d_dvc = NaN;
%             end 
%         end
%         
%         % conrad method.
%         function conrad(f,slope)
%             
%             if ~exist('slope','var')
%                 slope = 1;
%             end
% 
%             % Non-negative least-squares with relaxed OAR constraint
%             cvx_begin quiet
%                 variable x1(f.nBeamlts)
%                 minimize( sum_square(f.A_unif*x1 - f.d_unif) )
%                 subject to
%                     f.lb <= x1;
%                     for i = 1:f.nStructs
%                         for j = 1:f.structs{i}.nTerms
%                             if ~strcmp(f.structs{i}.terms{j}.type,'unif')
%                                 A = f.structs{i}.A;
%                                 d = f.structs{i}.terms{j}.d;
%                                 if strcmp(f.structs{i}.terms{j}.type,'udvc')
%                                     k = f.structs{i}.nVoxels - f.structs{i}.terms{j}.k;
%                                 else
%                                     k = f.structs{i}.terms{j}.k;
%                                 end
%                                 sum(pos(1 + slope*(A*x1 - d))) <= k;
%                             end
%                         end
%                     end
%             cvx_end
% 
%             % Non-negative least-squares with hard OAR constraint
%             cvx_begin quiet
%                 variable x2(f.nBeamlts)
%                 minimize( sum_square(f.A_unif*x2 - f.d_unif) )
%                 subject to
%                     f.lb <= x2;
%                     for i = 1:f.nStructs
%                         for j = 1:f.structs{i}.nTerms
%                             if ~strcmp(f.structs{i}.terms{j}.type,'unif')
%                                 [A,d] = calcConstMatSmall(f,x1,i,j);
%                                 A*x2 <= d;
%                             end
%                         end
%                     end
%             cvx_end
%             
%             f.x = x2;
%         end
%         
% 
%         
%         % Solve non-negative least-squares problem for x.
%         function projX(f)
%             
%             f.getd('full');
%             F = -f.A'*f.d;
%             fun = @(x)f.quadObj(x,f.H,F);
%             options.verbose = 0;
%             options.method = 'newton';
%             f.x = minConf_TMP(fun,f.x,f.lb,f.ub,options);
%         end
%         
%         
%         
%         % Plot objective function values.
%         function plotObj(f)
%             
%             figure()
%             myLines = lines;
%            
%             % Objective function
%             subplot(2,2,1)
%             plot(0:f.nIter,f.obj(1:f.nIter+1),'Color',[0.5,0.5,0.5])
%             xlabel('Iteration (k)')
%             ylabel('Objective Value')
%             
%             % Convergence of w variables
%             subplot(2,2,3)
%             plot(1:f.nIter,f.err(1:f.nIter),'Color',[0.5,0.5,0.5])
%             xlabel('Iteration (k)')
%             ylabel('Convergence Criteria')
%             
%             for i = 1:f.nStructs
%                 for j = 1:length(f.structs{i}.terms)
%                     % Objective function terms
%                     subplot(2,2,2), hold on
%                     plot(0:f.nIter,f.structs{i}.terms{j}.obj(1:f.nIter+1),'Color',myLines(i,:));
%             
%                     % Voxels under or over dose constraints
%                     if ~strcmp(f.structs{i}.terms{j}.type,'unif')
%                         subplot(2,2,4), hold on
%                         plot(0:f.nIter,f.structs{i}.terms{j}.vdiff(1:f.nIter+1),'Color',myLines(i,:));
%                         plot(0:f.nIter,f.structs{i}.terms{j}.wdiff(1:f.nIter+1),'--','Color',myLines(i,:));
%                     end
%                 end
%             end
%             subplot(2,2,2)
%             xlabel('Iteration (k)')
%             ylabel('Objective Terms')
%             legend(f.names,'Location','NorthEast')
%             
%             subplot(2,2,4)
%             xlabel('Iteration (k)')
%             ylabel('% Voxels Exceeding Dose')
%         end
%         
%         % Plot objective function values (fig 8).
%         function plotObjPaper(f)
%             
%             myLines = lines;
%            
%             % Objective function
%             figure(1)
%             subplot(3,1,1)
%             plot(0:f.nIter,f.obj(1:f.nIter+1),'Color',[0.5,0.5,0.5],'LineWidth',3)
%             f.adjustAxis(gca)
%             
%             % Objective terms
%             for i = 1:f.nStructs
%                 for j = 1:length(f.structs{i}.terms)
%                     figure(1)
%                     subplot(3,1,i+1)
%                     plot(0:f.nIter,f.structs{i}.terms{j}.obj(1:f.nIter+1),'Color',myLines(i,:),'LineWidth',3);
%                     f.adjustAxis(gca)
%             
%                     % Voxels under or over dose constraints
%                     if ~strcmp(f.structs{i}.terms{j}.type,'unif')
%                         figure(2), hold on
%                         subplot(2,1,2)
%                         plot(0:f.nIter,f.structs{i}.terms{j}.vdiff(1:f.nIter+1),'Color',myLines(i,:),'LineWidth',3);
%                         f.adjustAxis(gca)
%                         set(gca,'YTick',52:2:56);
%                     end
%                 end
%             end
%             
%             figure(2)
%             subplot(2,1,1)
%             plot(1:f.nIter,f.err(1:f.nIter),'Color',[0.5,0.5,0.5],'LineWidth',3)
%             f.adjustAxis(gca);
%         end   
%         
%         % Readjust axes limits.
%         function adjustAxis(~,g)
%             
%             axis tight
%             yVals = g.YLim;
%             yPad = 0.1*(yVals(2) - yVals(1));
%             g.YLim = [yVals(1)-yPad yVals(2)+yPad];
%             g.XTick = 0:50:200;
%             g.XTickLabels = {};
%             g.YTickLabels = {};
%             g.LineWidth = 2;      
%         end
%         
%         % Calculate and plot dose-volume histogram of solution.
%         function plotDVH(f)
%             
%             myLines = lines;
%             
%             % Calculate dose-volume histograms
%             doses = linspace(0,100,1000);
%             dvhInit = zeros(f.nStructs,length(doses));
%             dvhFinal = zeros(f.nStructs,length(doses));
%             for i = 1:f.nStructs
%                 doseInit = f.structs{i}.A*f.xInit;
%                 doseFinal = f.structs{i}.A*f.x;
%                 for j = 1:length(doses)
%                     dvhInit(i,j) = 100*sum(doseInit > doses(j))/f.structs{i}.nVoxels;
%                     dvhFinal(i,j) = 100*sum(doseFinal > doses(j))/f.structs{i}.nVoxels;
%                 end
%             end
%             
%             % Plot dose-volume histograms
%             figure(), hold on
%             
%             legendHandles = [];
%             legendNames = {};
%             for i = 1:f.nStructs
%                 for j = 1:length(f.structs{i}.terms)
%                     if ~strcmp(f.structs{i}.terms{j}.type,'unif') && f.structs{i}.terms{j}.percent == 0
%                         plot(f.structs{i}.terms{j}.dose,0,'p','MarkerFaceColor',[0.9290 0.6940 0.1250],...
%                             'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerSize',10);
%                     else
%                         if strcmp(f.structs{i}.terms{j}.type,'unif')
%                             percent = [0 100 100];
%                         elseif f.structs{i}.terms{j}.percent > 0
%                             percent = zeros(1,3);
%                             percent(2:3) = f.structs{i}.terms{j}.percent;
%                         end
%                         dose = zeros(1,3);
%                         dose(1:2) = f.structs{i}.terms{j}.dose;
%                         plot(dose,percent,':','Color',[0.4 0.4 0.4])
%                         plot(doses,dvhInit(i,:),'--','Color',myLines(i,:))
%                         if j == 1
%                             lineHandle = plot(doses,dvhFinal(i,:),'Color',myLines(i,:));
%                             lineName = f.structs{i}.name;
%                             legendHandles = [legendHandles lineHandle];
%                             legendNames = [legendNames, lineName];
%                         else
%                             plot(doses,dvhFinal(i,:),'Color',myLines(i,:))
%                         end
%                     end
%                 end
%                 
%                 % Annotations
%                 legend(legendHandles,legendNames)
%                 xlabel('Dose (Gy)')
%                 ylabel('Relative Volume (%)')
%                 ax = gca;
%                 ax.XLim = [0 doses(end)];
%                 ax.YLim = [0 100];
%                 box on
%                 axis square
%             end
%         end
%         
%         % Calculate and plot dose-volume histogram of solution (fig 9,11,12,13).
%         function plotDVHPaper(f)
%             
%             myLines = lines;
%             
%             % Calculate dose-volume histograms
%             doses = linspace(0,100,1000);
%             dvhInit = zeros(f.nStructs,length(doses));
%             dvhFinal = zeros(f.nStructs,length(doses));
%             for i = 1:f.nStructs
%                 doseInit = f.structs{i}.A*f.xInit;
%                 doseFinal = f.structs{i}.A*f.x;
%                 for j = 1:length(doses)
%                     dvhInit(i,j) = 100*sum(doseInit > doses(j))/f.structs{i}.nVoxels;
%                     dvhFinal(i,j) = 100*sum(doseFinal > doses(j))/f.structs{i}.nVoxels;
%                 end
%             end
%             
%             % Plot dose-volume histograms
%             for i = 1:f.nStructs
%                 figure(), hold on
%                 for j = 1:length(f.structs{i}.terms)
%                     if ~strcmp(f.structs{i}.terms{j}.type,'unif') && f.structs{i}.terms{j}.percent == 0
%                         plot(f.structs{i}.terms{j}.dose,0,'p','MarkerFaceColor',[0.9290 0.6940 0.1250],...
%                             'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerSize',10);
%                     else
%                         if strcmp(f.structs{i}.terms{j}.type,'unif')
%                             percent = [0 100 100];
%                         elseif f.structs{i}.terms{j}.percent > 0
%                             percent = zeros(1,3);
%                             percent(2:3) = f.structs{i}.terms{j}.percent;
%                         end
%                         dose = zeros(1,3);
%                         dose(1:2) = f.structs{i}.terms{j}.dose;
%                         plot(dose,percent,':','Color',[0.4 0.4 0.4],'LineWidth',3)
%                         plot(doses,dvhInit(i,:),'--','LineWidth',3,'Color',myLines(i,:))
%                         plot(doses,dvhFinal(i,:),'LineWidth',3,'Color',myLines(i,:))
%                     end
%                 end
%                 
%                 % Annotations
%                 ax = gca;
%                 ax.XLim = [0 doses(end)];
%                 ax.YLim = [0 100];
%                 ax.XTick = 0:20:100;
%                 ax.YTick = 0:20:100;
%                 ax.XTickLabel = {};
%                 ax.YTickLabel = {};
%                 ax.LineWidth = 2;
%                 box on
%                 axis square
%             end
%         end
%         
%         % Plot beamlet intensities.
%         function plotBeamlets(f)
%             
%             figure()
%             x = f.x;
%             
%             for i = 1:f.nAngles
%                 % Get x and y positions
%                 [linIdx,nx,ny] = f.getBeamlets(f.angles(i));
%                 C = zeros(nx,ny);
%                 C(linIdx) = 1;
%                 
%                 % Get beamlet intensities
%                 xTemp = x(1:length(linIdx));
%                 x = x(length(linIdx)+1:end);
%                 B = zeros(nx,ny);
%                 B(linIdx) = xTemp;
%                 B = B';
%                 
%                 % Plot beamlet intensities
%                 subplot(1,f.nAngles,i)
%                 imagesc(B), colormap gray
%                 title(sprintf('%d^\\circ',f.angles(i)),'Interpreter','tex')
%                 set(gca,'YDir','normal','XTick',[],'YTick',[])
%                 caxis([0 max(f.x)])
%                 axis square
%             end
%             p = get(subplot(1,f.nAngles,i),'Position');
%             cb = colorbar;
%             cb.Label.String = 'Beamlet Intensity (MU)';
%             set(subplot(1,f.nAngles,i),'Position',p);
%         end
%         
%         % Plot beamlet intensities (fig 10,14).
%         function plotBeamletsPaper(f)
%             
%             figure()
%             x = f.x;
%             
%             for i = 1:4
%                 % Get x and y positions
%                 [linIdx,nx,ny] = f.getBeamlets(f.angles(i));
%                 
%                 % Get beamlet intensities
%                 xTemp = x(1:length(linIdx));
%                 x = x(length(linIdx)+1:end);
%                 B = zeros(nx,ny);
%                 B(linIdx) = xTemp;
%                 B = B';
%                 
%                 % Plot beamlet intensities
%                 subplot(2,2,i)
%                 imagesc(B), colormap gray
%                 set(gca,'YDir','normal','XTick',[],'YTick',[])
%                 caxis([0 max(f.x)])
%                 axis square
%             end
%             
%             % Positioning
%             h = gcf;
%             pos = h.Position;
%             h.Position = [pos(1) pos(2) pos(3) pos(3)];
%             a = subplot(2,2,1);
%             a.Position = [0.1 0.5 0.3 0.3];
%             b = subplot(2,2,2);
%             b.Position = [0.45 0.5 0.3 0.3];
%             c = subplot(2,2,3);
%             c.Position = [0.1 0.15 0.3 0.3];
%             d = subplot(2,2,4);
%             d.Position = [0.45 0.15 0.3 0.3];
%             
%             % Colorbar
%             e = colorbar('southoutside','Ticks',0:1000:3000,'TickLabels',{},'LineWidth',2);
%             e.Position = [0.1    0.077    0.65    0.02700];
%         end
%         
%         % Get x and y positions for beamlets.
%         function [linIdx,nx,ny] = getBeamlets(~,angle)
%             
%             load(['Gantry' int2str(angle) '_Couch0_BEAMINFO.mat']);
%             xIdx = x - min(x) + 1;
%             yIdx = y - min(y) + 1;
%             nx = max(xIdx);
%             ny = max(yIdx);
%             linIdx = sub2ind([nx,ny],xIdx,yIdx);
%         end
%         
%         % Plot dose with slider for z position.
%         function plotDose(f)
%             
%             figure()
%             
%             % Get dose
%             Dose = reshape(f.D*f.x,184,184,90);
%             minDose = min(Dose(:));
%             maxDose = max(Dose(:));
%             str = 'z = %d';
%             c = [minDose maxDose];
%             
%             % Plot dose and contours
%             warning('off','MATLAB:contour:ConstantData');
%             hax = axes('Units','pixels');
%             imagesc(Dose(:,:,50),'AlphaData',Dose(:,:,50)~=0), hold on
%             for i = 1:length(f.mask)
%                 contour(f.mask{i}(:,:,50),1,'k');
%             end
%             
%             % Annotations
%             title(sprintf(str,50))
%             caxis(c);
%             cb = colorbar;
%             cb.Label.String = 'Dose (Gy)';
%             axis equal
%             axis off
%             hold off
%             
%             % Add slider
%             uicontrol('Style','slider',...
%                 'Min',1,'Max',90,'Value',50,...
%                 'Position',[200 20 120 20],...
%                 'Callback',{@f.whichSlice,hax,Dose,str,c});
%         end
%         
%         % Plot dose at slice 50 (fig 1,5,10,14).
%         function plotDosePaper(f)
%             
%             figure()
%             idx1 = 40:126;
%             idx2 = 23:152;
%             
%             % Get CT slice
%             ct = dicomread('CT.2.16.840.1.113662.2.12.0.3173.1271873797.276');
%             ct = double(imresize(ct,[184,184]));
%             ct50 = ct(idx1,idx2);
%             ctShift = ct50 - min(ct50(:));
%             ctShiftScale = ctShift/max(ctShift(:));
%             CT50 = repmat(ctShiftScale,[1 1 3]);
%             
%             % Get Dose
%             Dose = reshape(f.D*f.x,184,184,90);
%             Dose50 = Dose(idx1,idx2,50);
%             
%             % Plot CT
%             body50 = f.mask{end}(idx1,idx2,50);
%             imagesc(CT50), hold on
%             
%             % Plot dose
%             imagesc(Dose50,'AlphaData',0.3*(body50~=0))
%             contour(Dose50,0:10:100,'LineWidth',2);
%             
%             % Plot organ contours
%             for i = 1:length(f.mask)-1
%                contour(f.mask{i}(idx1,idx2,50),1,'k','LineWidth',2); 
%             end
%             
%             % Annotations
%             caxis([min(Dose50(:)) max(Dose50(:))]);
%             axis equal
%             axis off
%             
%             % Colorbar
%             colorbar('southoutside','Ticks',0:20:100,'TickLabels',{},'LineWidth',2)
%         end
%         
%         % Callback function for plotDose() slider.
%         function whichSlice(f,hObj,~,~,Dose,str,c)
%             
%             % Plot dose and contours
%             z = round(get(hObj,'Value'));
%             imagesc(Dose(:,:,z),'AlphaData',Dose(:,:,z)~=0), hold on
%             for i = 1:length(f.mask)
%                 contour(f.mask{i}(:,:,z),1,'k');
%             end
%             
%             % Annotations
%             title(sprintf(str,z))
%             caxis(c);
%             cb = colorbar;
%             cb.Label.String = 'Dose (Gy)';
%             axis equal
%             axis off   
%             hold off
%         end
%         
%         % Save results
%         function saveResults(f,str)
%             
%             % Input parameters
%             pars.structs = f.structs;
%             pars.angles = f.angles;
%             pars.lambda = f.lambda;
%             pars.maxIter = f.maxIter;
%             pars.xInit = f.xInit;
%             pars.overlap = f.overlap;
%             pars.tol = f.tol;
%             
%             % Solution variables
%             pars.x = f.x;
%             pars.obj = f.obj;
%             pars.nIter = f.nIter;
%             save(str,'pars');
%         end
        end
    end
end

function D = getD(angles)
% GETD Get full beamlet-to-voxel matrix.

    temp = [];
    for ii = angles
        load(['Gantry' int2str(ii) '_Couch0_D.mat']);
        temp = [temp D];
    end
    D = temp;
end

function structs = getStructVars(structs,overlap,D)
% GETSTRUCTVARS Get structure-specific variables.

    vPrev = [];
    for ii = 1:length(structs)
        load([structs{ii}.name '_VOILIST.mat']);
        if ~overlap
            [v,vPrev] = removeOverlap(v,vPrev); 
        end
        structs{ii}.A = D(v,:);
        structs{ii}.terms = getTermVars(structs{ii}.terms,length(v));
    end
end

function [v,vPrev] = removeOverlap(v,vPrev)
% REMOVEOVERLAP Remove overlapping voxels from body structure.

   if isempty(vPrev)
       vPrev = v;
   else
       v = setdiff(v,vPrev);
       vPrev = union(v,vPrev);
   end 
end

function terms = getTermVars(terms,nVoxels)
% GETTERMVARS Get term-specific variables.

    for ii = 1:length(terms)
        terms{ii}.d = terms{ii}.dose*ones(nVoxels,1);
        terms{ii}.step = nVoxels/terms{ii}.weight;
        if ~strcmp(terms{ii}.type,'unif')
            % number of voxels allowed to be < or > dose
            terms{ii}.k = floor(terms{ii}.percent*nVoxels/100);
        end
    end
end

function names = getNames(structs)
% GETNAMES Get body structure names.

    nStructs = length(structs);
    names = cell(1,nStructs);
    for ii = 1:nStructs
        names{ii} = structs{ii}.name;
    end
end

function mask = getMaskStruct(names,overlap)
% GETMASKSTRUCT Get body structure contours for all organs.

    vPrev = [];
    for ii = 1:length(names)
       load([names{ii} '_VOILIST.mat']);
       if ~overlap
           [v,vPrev] = removeOverlap(v,vPrev); 
       end
       mask{ii} = getMask(v);
    end
    if ~any(strcmp(names,'BODY'))
        load('BODY_VOILIST.mat');
        mask{ii+1} = getMask(v);
    end
end

function mask = getMask(v)
% GETMASK Get body structure contour for one organ.

    mask = zeros(184*184*90,1);
    mask(v) = 1;
    mask = reshape(mask,184,184,90);
end

function [A,lb,ub] = getA(structs,lambda,varargin)
% GETA Get stacked beamlet-to-voxel matrix and beamlet bounds.
%
%   [A,H,lb,ub] = getA(structList,lambda)
%       Get output for all structures and terms.
%
%   [A,H,lb,ub] = getA(structList,lambda,type)
%       Get output for structures and terms specified by type:
%           'full' - All structures and terms
%           'unif' - Only structures and terms with uniform targets

    % Get type of matrix
    if nargin > 2
        type = varargin{1};
    else
        type = 'full';
    end
    matFull = strcmp(type,'full');
    matUnif = strcmp(type,'unif');
            
    % Add terms
    A = [];
    for ii = 1:length(structs)
        for jj = 1:length(structs{ii}.terms)
            termUnif = strcmp(structs{ii}.terms{jj}.type,'unif');
            if matFull || (matUnif && termUnif)
                weight = structs{ii}.terms{jj}.weight;
                nVoxels = size(structs{ii}.A,1);
                temp = sqrt(weight/nVoxels)*structs{ii}.A;
                A = [A; temp];
            end
        end
    end
    
    % Add regularization
    if lambda > 0
        A = [A; sqrt(lambda)*eye(size(A,2))];
    end
    
    % Create beamlet bounds
    lb = zeros(size(A,2),1);
    ub = inf(size(lb));
end

function d = getd(structs,lambda,varargin)
% GETD Get stacked dose vector.
%
%   d = getd(structs)
%       Get output for all structures and terms.
%
%   d = getd(structs,type)
%       Get output for structurs and terms specified by type:
%           'full' - All structures and terms
%           'unif' - Only structurs and terms with uniform targets

    % Get type of vector
    if nargin > 2
        type = varargin{1};
    else
        type = 'full';
    end
    vecFull = strcmp(type,'full');
    vecUnif = strcmp(type,'unif');
    
    % Add terms
    d = [];
    for ii = 1:length(structs)
        nVoxels = size(structs{ii}.A,1);
        for jj = 1:length(structs{ii}.terms)
            termUnif = strcmp(structs{ii}.terms{jj}.type,'unif');
            if vecFull || (vecUnif && termUnif)
                weight = structs{ii}.terms{jj}.weight;
                temp = sqrt(weight/nVoxels)*structs{ii}.terms{jj}.d;
                if ~termUnif
                    w = structs{ii}.terms{jj}.w;
                    temp = temp + sqrt(weight/nVoxels)*w;
                end
                d = [d; temp];
            end
        end
    end
    
    % Add regularization
    if lambda > 0
        d = [d; zeros(size(structs{1}.A,2),1)];
    end
end

function x = getInitX(structs,lambda)
% GETINITX Initialze beamlets.

    [A,lb,ub] = getA(structs,lambda,'unif');
    d = getd(structs,lambda,'unif');
    options = optimoptions(@lsqlin,'Display','off');
    x = lsqlin(A,d,[],[],[],[],lb,ub,[],options);
end

% need to test...
function w = projW(w,k)
% PROJW Project w onto the set satisfying ||max(0,w)||_0 <= k.

    idxPos = w > 0;
    if sum(idxPos) > k
        wPos = w(idxPos);
        [~,idxSort] = sort(wPos,'descend');
        wPos(idxSort(k+1:end)) = 0;
        w(idxPos) = wPos;
    end
end
