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
    %       weight: coefficient of the term in the objective function
    %       dose: dose in Gy
    %       percent: include if type is 'ldvc' or 'udvc':
    %           * 'ldvc': No more than p% receives less than d Gy
    %           * 'udvc': No more than p% receives more than d Gy
    %
    %   Written to work with the CORT prostate tumor dataset, but could be
    %   modified to work with other datasets. 
    
    properties (SetAccess = private)
        structs            % Body structures
        angles = 0:52:358; % Gantry angles
        lambda = 1e-8;     % L2 regularization coefficient
        overlap = false;   % Allow overlaps in structures
        nStructs           % Number of body structures
        nAngles            % Number of angles
        nBeamlts           % Number of beamlets
    end
    
    properties (Access = private)
        D     % Full beamlet-to-voxel matrix
        A     % Stacked beamlet-to-voxel matrix
        Au    % Stacked beamlet-to-voxel matrix for uniform target terms
        du    % Stacked dose vector for uniform target terms
        lb    % Lower bound for beamlets
        ub    % Upper bound for beamlets
        names % Body structure names
        mask  % Body structure contours for plotting     
    end

    properties
        x0             % Initial beamlet intensities
        x              % Final beamlet intensities
        obj            % Objective function values
        wDiff          % Convergence criteria
        nIter          % Number of iterations used
        tol = 1e-3;    % Stopping tolerance
        maxIter = 500; % Maximum number of iterations
    end
    
    methods
        function prob = FluenceMapOpt(structs,varargin)
            % FLUENCEMAPOPT Initialize problem.
            %
            %   prob = FluenceMapOpt(structs)
            %       Initialize problem with default parameters.
            %
            %   prob = FluenceMapOpt(structs,xInit,angles,lambda,overlap,tol,maxIter)
            %       Initialize problem with optional arguments.
        
            % Set input variables
            if nargin == 0
                error('Not enough input arguments.')
            end
            if ~iscell(structs)
                error('Invalid input for `structs`.')
            end
            prob.setInputVars(varargin,nargin-1);
            prob.nStructs = length(structs);
            prob.nAngles = length(prob.angles);
            
            % Comput internal variables
            [prob.D,prob.nBeamlts] = FluenceMapOpt.getD(prob.angles);
            prob.structs = FluenceMapOpt.getStructVars(structs,...
                prob.nStructs,prob.overlap,prob.D);
            prob.names = FluenceMapOpt.getNames(prob.structs,prob.nStructs);
            prob.mask = FluenceMapOpt.getMaskStruct(prob.names,prob.overlap);
            [prob.A,prob.lb,prob.ub] = FluenceMapOpt.getA(prob.structs,...
                prob.lambda,prob.nStructs,prob.nBeamlts);
            [prob.Au,~,~] = FluenceMapOpt.getA(prob.structs,prob.lambda,...
                prob.nStructs,prob.nBeamlts,'unif');
            prob.du = FluenceMapOpt.getd(prob.structs,prob.lambda,...
                prob.nStructs,prob.nBeamlts,'unif');
            
            % Compute initial beamlets
            if nargin > 1
                prob.x0 = varargin{1};
            else
                prob.x0 = FluenceMapOpt.projX(prob.Au,prob.du,prob.lb,prob.ub);
            end
        end
        
        % need to test...
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
            for kk = 1:prob.maxIter
                
                % Update x and w vectors
                prob.updateX();
                wDiffSum = 0;
                for ii = 1:prob.nStructs
                    for jj = 1:prob.structs{ii}.nTerms
                        if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                            wDiffSum = wDiffSum + prob.updateW(ii,jj);
                        end
                    end
                end
                
                % Calculate objective
                prob.nIter = kk;
                prob.calcObj(kk,print);
                prob.wDiff(kk) = wDiffSum;
                if print
                    fprintf(', wDiff: %7.4e\n',wDiffSum);
                end
                
                % Check convergence
                if wDiffSum <= prob.tol
                    break
                end
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
    
    methods (Hidden)
        function setInputVars(prob,args,nArgs)
            % SETINPUTVARS Set input variables.
            varNames = {'angles','lambda','overlap','tol','maxIter'};
            for ii = 1:length(varNames)
               if nArgs > ii
                   if ~isempty(args{ii+1})
                       prob.(varNames{ii}) = args{ii+1};
                   end
               end
            end
        end 
        
        function initProb(prob,print)
            % INITPROB Initialize x, w, and objective values.
            prob.x = prob.x0;
            prob.initW();
            prob.initObj();
            prob.calcObj(0,print);
        end

        function initW(prob)
            % INITW Initialize w vectors for dose-volume constraint terms.
            for ii = 1:prob.nStructs
                for jj = 1:prob.structs{ii}.nTerms
                    if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        initDose = prob.structs{ii}.A*prob.x0;
                        res = initDose - prob.structs{ii}.terms{jj}.d;
                        s = strcmp(prob.structs{ii}.terms{jj}.type,'ldvc');
                        k = prob.structs{ii}.terms{jj}.k;
                        w = FluenceMapOpt.projW((-1)^s*res,k);
                        prob.structs{ii}.terms{jj}.w = w;
                   end
                end
            end
        end
        
        function initObj(prob)
            % INITOBJ Initialize objective function values
            myZeros = zeros(1,prob.maxIter+1);
            prob.obj = myZeros;
            prob.wDiff = myZeros(1:end-1);
            for ii = 1:prob.nStructs
                for jj = prob.structs{ii}.nTerms
                    prob.structs{ii}.terms{jj}.obj = myZeros;
                    if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        prob.structs{ii}.terms{jj}.resPos = myZeros;
                        prob.structs{ii}.terms{jj}.wPos = myZeros;
                    end
                end
            end
        end
        
        function calcObj(prob,iter,print)
            % CALCOBJ Calculate and print objective function values.
            for ii = 1:prob.nStructs
                nVoxels = prob.structs{ii}.nVoxels;
                for jj = 1:prob.structs{ii}.nTerms
                    weight = prob.structs{ii}.terms{jj}.weight;
                    dose = prob.structs{ii}.A*prob.x;
                    res = dose - prob.structs{ii}.terms{jj}.d;
                    if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        s = strcmp(prob.structs{ii}.terms{jj}.type,'ldvc');
                        resPos = 100*sum((-1)^s*res > 0)/nVoxels;
                        prob.structs{ii}.terms{jj}.dPos(iter+1) = resPos;
                        wPos = 100*sum(prob.structs{ii}.terms{jj}.w > 0)/nVoxels;
                        prob.structs{ii}.terms{jj}.wPos(iter+1) = wPos;
                        res = res - prob.structs{ii}.terms{jj}.w; 
                    end
                    termObj = weight*norm(res)^2/(2*nVoxels);
                    prob.structs{ii}.terms{jj}.obj(iter+1) = termObj;
                    prob.obj(iter+1) = prob.obj(iter+1) + termObj;
                end
            end
            prob.obj(iter+1) = prob.obj(iter+1) + prob.lambda*norm(prob.x)^2/2;
            if print
                fprintf('iter: %d, obj: %7.4e',iter,prob.obj(iter+1));
                if iter == 0
                    fprintf('\n');
                end
            end 
        end
        
        function updateX(prob)
            % UPDATEX Update beamlet intensities.
            d = FluenceMapOpt.getd(prob.structs,prob.lambda,...
                prob.nStructs,prob.nBeamlts);
            prob.x = FluenceMapOpt.projX(prob.A,d,prob.lb,prob.ub);
        end
        
        % need to test...
        function wDiff = updateW(prob,ii,jj)
            % UPDATEW Update proxy variable.
            
            % Grab variables
            s = strcmp(prob.structs{ii}.terms{jj}.type,'ldvc');
            k = prob.structs{ii}.terms{jj}.k;
            step = prob.structs{ii}.terms{jj}.step;
            coeff = step/prob.structs{ii}.nVoxels;
            
            % Calculate gradient step
            dose = prob.structs{ii}.A*prob.x;
            res = (-1)^s*(dose - prob.structs{ii}.terms{jj}.d);
            wPrev = prob.structs{ii}.terms{jj}.w;
            wStep = wPrev + coeff*(res - wPrev);
            
            % Project onto set ||(w)_+||_0 <= k
            wProj = FluenceMapOpt.projW(wStep,k);
            wDiff = norm(wProj - wPrev)/step;
            prob.structs{ii}.terms{jj}.w = wProj;
        end
    end
    
    methods (Hidden, Static)
        function [D,nBeamlts] = getD(angles)
            % GETD Get full beamlet-to-voxel matrix.
            temp = [];
            for ii = angles
                load(['Gantry' int2str(ii) '_Couch0_D.mat']);
                temp = [temp D];
            end
            D = temp;
            nBeamlts = size(D,2);
        end

        function structs = getStructVars(structs,nStructs,overlap,D)
            % GETSTRUCTVARS Get structure-specific variables.
            vPrev = [];
            for ii = 1:nStructs
                load([structs{ii}.name '_VOILIST.mat']);
                if ~overlap
                    [v,vPrev] = FluenceMapOpt.removeOverlap(v,vPrev); 
                end
                structs{ii}.A = D(v,:);
                structs{ii}.nVoxels = length(v);
                structs{ii}.nTerms = length(structs{ii}.terms);
                structs{ii}.terms = FluenceMapOpt.getTermVars(structs{ii}.terms,...
                    structs{ii}.nTerms,structs{ii}.nVoxels);
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

        function terms = getTermVars(terms,nTerms,nVoxels)
            % GETTERMVARS Get term-specific variables.
            for ii = 1:nTerms
                terms{ii}.d = terms{ii}.dose*ones(nVoxels,1);
                terms{ii}.step = nVoxels/terms{ii}.weight;
                if ~strcmp(terms{ii}.type,'unif')
                    % number of voxels allowed to be < or > dose
                    terms{ii}.k = floor(terms{ii}.percent*nVoxels/100);
                end
            end
        end

        function names = getNames(structs,nStructs)
            % GETNAMES Get body structure names.
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
                   [v,vPrev] = FluenceMapOpt.removeOverlap(v,vPrev); 
               end
               mask{ii} = FluenceMapOpt.getMask(v);
            end
            if ~any(strcmp(names,'BODY'))
                load('BODY_VOILIST.mat');
                mask{ii+1} = FluenceMapOpt.getMask(v);
            end
        end

        function mask = getMask(v)
            % GETMASK Get body structure contour for one organ.
            mask = zeros(184*184*90,1);
            mask(v) = 1;
            mask = reshape(mask,184,184,90);
        end

        function [A,lb,ub] = getA(structs,lambda,nStructs,nBeamlts,varargin)
            % GETA Get stacked beamlet-to-voxel matrix and beamlet bounds.
            %
            %   [A,H,lb,ub] = getA(structs,lambda,nStructs,nBeamlts)
            %       Get output for all structures and terms.
            %
            %   [A,H,lb,ub] = getA(structs,lambda,nStructs,nBeamlts,type)
            %       Get output for structures and terms specified by type:
            %           * 'full': All structures and terms
            %           * 'unif': Structures and terms with uniform targets
            if nargin > 4
                type = varargin{1};
            else
                type = 'full';
            end
            matFull = strcmp(type,'full');
            matUnif = strcmp(type,'unif');

            % Add terms
            A = [];
            for ii = 1:nStructs
                for jj = 1:structs{ii}.nTerms
                    termUnif = strcmp(structs{ii}.terms{jj}.type,'unif');
                    if matFull || (matUnif && termUnif)
                        weight = structs{ii}.terms{jj}.weight;
                        nVoxels = structs{ii}.nVoxels;
                        temp = sqrt(weight/nVoxels)*structs{ii}.A;
                        A = [A; temp];
                    end
                end
            end

            % Add regularization
            if lambda > 0
                A = [A; sqrt(lambda)*eye(nBeamlts)];
            end

            % Create beamlet bounds
            lb = zeros(nBeamlts,1);
            ub = inf(nBeamlts,1);
        end

        function d = getd(structs,lambda,nStructs,nBeamlts,varargin)
            % GETD Get stacked dose vector.
            %
            %   d = getd(structs,lambda,nStructs,nBeamlts)
            %       Get output for all structures and terms.
            %
            %   d = getd(structs,lambda,nStructs,nBeamlts,type)
            %       Get output for structurs and terms specified by type:
            %           * 'full':  All structures and terms
            %           * 'unif': Structures and terms with uniform targets

            % Get type of vector
            if nargin > 4
                type = varargin{1};
            else
                type = 'full';
            end
            vecFull = strcmp(type,'full');
            vecUnif = strcmp(type,'unif');

            % Add terms
            d = [];
            for ii = 1:nStructs
                nVoxels = structs{ii}.nVoxels;
                for jj = 1:structs{ii}.nTerms
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
                d = [d; zeros(nBeamlts,1)];
            end
        end

        function x = projX(A,d,lb,ub)
            % PROJX Solve non-negative least-squares problem for beamlets.
            options = optimoptions(@lsqlin,'Display','off');
            x = lsqlin(A,d,[],[],[],[],lb,ub,[],options);
        end

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
    end
end
