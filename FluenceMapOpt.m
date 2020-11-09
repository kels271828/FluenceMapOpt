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
        structs               % Body structures
        angles = 0:52:358;    % Gantry angles
        overlap = false;      % Allow overlaps in structures
        lambda = 1e-8;        % L2 regularization coefficient
        nnls = 'minConf_TMP'; % Method to compute NNLS problem
        nStructs              % Number of body structures
        nDVC                  % Number of dose-volume constraints
        nAngles               % Number of angles
        nBeamlets             % Number of beamlets
        obj                   % Objective function values
        wDiff                 % Convergence criteria
        nIter                 % Number of iterations used
        time                  % Time to compute solution (seconds)
    end
    
    properties (Access = private)
        D     % Full beamlet-to-voxel matrix
        A     % Stacked beamlet-to-voxel matrix
        H     % Stacked beamlet-to-voxel Hessian
        Au    % Stacked beamlet-to-voxel matrix for uniform target terms
        Hu    % Stacked beamlet to voxel Hessian for uniform target terms
        du    % Stacked dose vector for uniform target terms
        As    % Stacked beamlet-to-voxel matrix for model with slack
        Hs    % Stacked beamlet-to-voxel Hessian for model with slack
        lb    % Lower bound for beamlets
        ub    % Upper bound for beamlets
        lbs   % Lower bound for beamlets and slack
        ubs   % Upper bound for beamlets and slack
        names % Body structure names
        mask  % Body structure contours for plotting     
    end

    properties
        x0              % Initial beamlet intensities
        x               % Final beamlet intensities
        tol = 1e-3;     % Stopping tolerance
        maxIter = 1000; % Maximum number of iterations
    end
    
    methods
        function prob = FluenceMapOpt(structs,varargin)
            % FLUENCEMAPOPT Initialize problem.
            %
            %   prob = FluenceMapOpt(structs)
            %       Initialize problem with default parameters.
            %
            %   prob = FluenceMapOpt(structs,ProbSpec)
            %       Initialize problem with optional arguments.
            %
            %   Example:
            %       prob = FluenceMapOpt(structs,...
            %           'angles',0:52:358,...
            %           'overlap',false,...
            %           'lambda',1e-8,...
            %           'nnls','minConf_TMP',...
            %           'x0',zeros(986,1),...
            %           'tol',1e-3,...
            %           'maxIter',1000);
        
            % Set input variables
            if nargin == 0
                error('Not enough input arguments.')
            end
            if ~iscell(structs)
                error('Invalid input for `structs`.')
            end
            prob.setInputVars(varargin);
            prob.nStructs = length(structs);
            prob.nAngles = length(prob.angles);
            
            % Comput internal variables
            [prob.D,prob.nBeamlets] = FluenceMapOpt.getD(prob.angles);
            [prob.structs,prob.nDVC] = FluenceMapOpt.getStructVars(structs,...
                prob.nStructs,prob.overlap,prob.D);
            prob.names = FluenceMapOpt.getNames(prob.structs,prob.nStructs);
            prob.mask = FluenceMapOpt.getMaskStruct(prob.names,prob.overlap);
            [prob.A,prob.H,prob.lb,prob.ub] = prob.getA('full');
            [prob.Au,prob.Hu,~,~] = prob.getA('unif');
            prob.du = prob.getd('unif');
            
            % Compute initial beamlets
            if isempty(prob.x0)
                prob.x0 = prob.projX('unif');
            end
            prob.x = prob.x0;
        end
        
        function updateStructs(prob,structs,x0)
            % UPDATESTRUCTS Update structure weights, doses, or percents.
            if nargin == 1
                error('Not enough input arguments.')
            end
            if ~iscell(structs)
                error('Invalid input for `structs`.')
            end
            for ii = 1:prob.nStructs
                prob.structs{ii}.terms = FluenceMapOpt.getTermVars(...
                    structs{ii}.terms,prob.structs{ii}.nTerms,...
                    prob.structs{ii}.nVoxels);
            end
            [prob.A,prob.H,~,~] = prob.getA('full');
            [prob.Au,prob.Hu,~,~] = prob.getA('unif');
            prob.du = prob.getd('unif');
            if nargin < 3
                prob.x0 = prob.projX('unif');
            else
                prob.x0 = x0;
            end
            prob.x = prob.x0;
        end
        
        function calcBeams(prob,print)
            % CALCBEAMLETS Calculate beamlet intensities.
            if nargin == 1
                print = true;
            end
            
            % Fluence map optimization
            tic;
            prob.initProb(print);
            for kk = 1:prob.maxIter
                
                % Update x and w vectors
                prob.x = prob.projX('full');
                wDiffSum = 0;
                for ii = 1:prob.nStructs
                    for jj = 1:prob.structs{ii}.nTerms
                        if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                            wDiffSum = wDiffSum + prob.updateW(ii,jj);
                        end
                    end
                end
                wDiffSum = wDiffSum/prob.nDVC;
                
                % Calculate objective
                prob.nIter = kk;
                prob.calcObj(kk,print);
                prob.wDiff(kk) = wDiffSum;
                if print
                    fprintf(', wDiff: %7.4e\n',wDiffSum);
                end
                
                % Check convergence
                if wDiffSum <= prob.tol
                    prob.obj = prob.obj(1:prob.nIter+1);
                    prob.wDiff = prob.wDiff(1:prob.nIter);
                    break
                end
            end
            prob.x = prob.projX('full');
            prob.time = toc;
        end
        
        function calcBeamsConvex(prob,print)
            % CONVRELAX Approach inspired by Fu paper.
            if nargin == 1
                print = true;
            end
            
            % Get number of dose-volume constraint terms
            numA = 0;
            for ii = 1:prob.nStructs
                for jj = 1:prob.structs{ii}.nTerms
                    if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        numA = numA + 1;
                    end
                end
            end
           
            % Fluence map optimization
            tic;
            if print
                cvx_begin
            else
                cvx_begin quiet
            end
            variables xRelax(prob.nBeamlets) a(numA)
            minimize(sum_square(prob.Au*xRelax - prob.du))
            subject to
                prob.lb <= xRelax;
                zeros(numA,1) <= a;
                idxA = 1;
                for ii = 1:prob.nStructs
                    At = prob.structs{ii}.A;
                    nVoxels = prob.structs{ii}.nVoxels;
                    for jj = 1:prob.structs{ii}.nTerms
                        if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                            d = prob.structs{ii}.terms{jj}.d;
                            k = prob.structs{ii}.terms{jj}.k;
                            s = strcmp(prob.structs{ii}.terms{jj}.type,'ldvc');
                            res = (At*xRelax - d);
                            termSum = sum(pos(a(idxA) + (-1)^s*res));
                            termSum <= a(idxA)*k;
                            idxA = idxA + 1;
                        end
                    end
                end
            cvx_end
            prob.x = xRelax;
            prob.time = toc;
        end
        
        function calcBeamsIter(prob,print)
            % CALCBEAMSITER Approach inspired by Llacer paper.
            if nargin == 1
                print = true;
            end
            
            % Fluence map optimization
            tic;
            prob.x = prob.x0;
            for kk = 1:prob.maxIter
                xOld = prob.x;
                grad = prob.getIterGrad(prob.x);
                prob.x = prob.x - grad;
                prob.x(prob.x < 0) = 0;
                xDiff = norm(xOld - prob.x)/prob.nBeamlets;
                
                % Calculate objective
                prob.nIter = kk;
                if print
                    obj = prob.getIterObj(prob.x);        
                    fprintf('iter: %d, obj: %7.4e, xDiff: %7.4e\n',...
                        kk,obj,xDiff);
                end
                
                % Check convergence
                if xDiff <= prob.tol
                    break
                end
            end
            prob.time = toc;
        end
        
        function calcBeamsSlack(prob,print)
            % CALCBEAMSLACK Approach inspired by Zhang paper.
            if nargin == 1
                print = true;
            end
            if isempty(prob.As)
                [prob.As,prob.Hs,prob.lbs,prob.ubs] = prob.getA('slack');
            end
            
            % Fluence map optimization
            tic;
            prob.initU();
            y = zeros(size(prob.As,2),1);
            for kk = 1:prob.maxIter
                % Update x and u vectors
                y = prob.projX('slack',y);
                prob.x = y(1:prob.nBeamlets);
                uDiffSum = 0;
                for ii = 1:prob.nStructs
                    for jj = 1:prob.structs{ii}.nTerms
                        if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                            uDiffSum = uDiffSum + prob.updateU(ii,jj);
                        end
                    end
                end
                uDiffSum = uDiffSum/prob.nDVC;
                
                % Calculate objective
                prob.nIter = kk;
                if print
                    ds = prob.getd('slack');
                    obj = norm(prob.As*y - ds)^2/2;
                    fprintf('iter: %d, obj: %7.4e, uDiff: %7.4e\n',...
                        kk,obj,uDiffSum);
                end

                % Check convergence
                if uDiffSum <= prob.tol
                    break
                end
            end
            y = prob.projX('slack',y);
            prob.x = y(1:prob.nBeamlets);
            prob.time = toc;
        end
        
        function calcBeamsPolish(prob,x,print)
            % CALCBEAMSPOLISH Approach inspired by Saberian paper.
            %
            %   Display options: 
            %       * 'off', 'none', 'final', 'iter'
            %       * 'iter-detailed', 'final-detailed'
            if nargin == 2
                print = 'iter';
            else
                if ~print
                    print = 'off';
                end
            end
            tic;
            f = -prob.Au'*prob.du;
            [Ac,dc] = prob.getConstraints(x);
            options = optimoptions(@quadprog,'Display',print);
            prob.x = quadprog(prob.Hu,f,Ac,dc,[],[],prob.lb,prob.ub,[],options);
            prob.time = toc;
        end
        
        function plotObj(prob)
            % PLOTOBJ Plot objective function values.
            
            % Objective function
            figure()
            subplot(1,2,1)
            plot(0:prob.nIter,prob.obj,'LineWidth',2)
            xlabel('Iteration (k)')
            ylabel('Objective Value')
            
            % Convergence of proxy variables
            subplot(1,2,2)
            plot(1:prob.nIter,prob.wDiff,'LineWidth',2)
            xlabel('Iteration (k)')
            ylabel('Convergence Criteria')
        end
        
        function plotDVH(prob,legendNames)
            % PLOTDVH Plot dose-volume histograms for initial and final dose.
            if nargin == 1
                legendNames = prob.names;
            end
            
            % Compute curves and initialize
            [doses,dvhInit] = prob.calcDVH(prob.x0);
            [~,dvhFinal] = prob.calcDVH(prob.x);
            myLines = lines;
            legendHandles = [];
            figure(), hold on
            
            % Plot dose-volume histograms
            for ii = 1:prob.nStructs
                for jj = 1:prob.structs{ii}.nTerms
                    % Plot targets/constraints
                    prob.plotConstraints(ii,jj);
                end
                % Plot dvh curves
                plot(doses,dvhInit(ii,:),'--','Color',myLines(ii,:),...
                    'LineWidth',2)
                dvhHandle = plot(doses,dvhFinal(ii,:),...
                    'Color',myLines(ii,:),'LineWidth',2);
                legendHandles = [legendHandles dvhHandle];
            end
                
            % Annotations
            legend(legendHandles,legendNames,'Location','northeastoutside')
            xlabel('Dose (Gy)')
            ylabel('Relative Volume (%)')
            ax = gca;
            ax.XLim = [0 doses(end)];
            ax.YLim = [0 100];
            box on
            axis square
        end
        
        function plotDVHPaper(prob,xMat,unif)
            % PLOTDVHPAPER Plot dose-volume histogram of solution.
            if nargin == 2
                unif = true;
            end
            nX = size(xMat,2);
            
            % Compute curves and initialize
            dvhMat = [];
            for ii = 1:nX
                [doses,dvh] = prob.calcDVH(xMat(:,ii));
                dvhMat = cat(3,dvhMat,dvh);
            end
            myColors = lines;
            if nX == 2
                myLines = {'--','-'}';
            elseif nX == 3
                myLines = {'--','-',':'};
            else
                myLines = {'--','-',':','-.'};
            end            
            
            % Plot dose-volume histograms
            for ii = 1:prob.nStructs
                figure(), hold on
                % Plot targets/constraints
                for jj = 1:length(prob.structs{ii}.terms)
                    if unif || ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        prob.plotConstraints(ii,jj);
                    end
                end
                % Plot dvh curves
                for kk = 1:nX
                    if jj == prob.structs{ii}.nTerms
                        plot(doses,dvhMat(ii,:,kk),myLines{kk},...
                            'Color',myColors(1,:),'LineWidth',2);
                    end        
                end
                
                % Annotations
                ax = gca;
                ax.XLim = [0 doses(end)];
                ax.YLim = [0 100];
                ax.XTick = 0:20:100;
                ax.YTick = 0:20:100;
                ax.XTickLabel = {};
                ax.YTickLabel = {};
                ax.LineWidth = 2;
                box on
                axis square
            end
        end
        
        function compareDVH(prob,xMat,legendNames)
            % COMPAREDVH Plot dose-volume histograms for multiple solutions.
            nX = size(xMat,2);
            if nargin == 2
                legendNames = cell(1,nX);
                for ii = 1:nX
                    legendNames{ii} = sprintf('x%d',ii);
                end
            end                    
            
            % Compute curves and initialize
            dvhMat = [];
            for ii = 1:nX
                [doses,dvh] = prob.calcDVH(xMat(:,ii));
                dvhMat = cat(3,dvhMat,dvh);
            end
            myLines = lines;
            legendHandles = [];
            figure(), hold on
            
            % Plot dose-volume histograms
            for ii = 1:prob.nStructs
                for jj = 1:prob.structs{ii}.nTerms
                    % Plot targets/constraints
                    prob.plotConstraints(ii,jj);
                   
                    % Plot dvh curves
                    for kk = 1:nX
                        if ii == 1 && jj == 1
                            dvhHandle = plot(doses,dvhMat(ii,:,kk),...
                                'Color',myLines(kk,:),'LineWidth',2);
                            legendHandles = [legendHandles dvhHandle];
                        else
                            plot(doses,dvhMat(ii,:,kk),...
                                'Color',myLines(kk,:),'LineWidth',2)
                        end 
                    end
                end
            end
                
            % Annotations
            legend(legendHandles,legendNames,'Location','northeastoutside')
            xlabel('Dose (Gy)')
            ylabel('Relative Volume (%)')
            ax = gca;
            ax.XLim = [0 doses(end)];
            ax.YLim = [0 100];
            box on
            axis square
        end
        
        function plotBeams(prob)
            %PLOTBEAMS Plot beamlet intensities.
            figure()
            xRemain = prob.x;
            for ii = 1:prob.nAngles
                % Get beamlet intensities
                [idx,nX,nY] = FluenceMapOpt.getBeamCoords(prob.angles(ii));
                xCurrent = xRemain(1:length(idx));
                xRemain = xRemain(length(idx)+1:end);
                beam = zeros(nX,nY);
                beam(idx) = xCurrent;
                beam = beam';
                
                % Plot beamlet intensities
                subplot(1,prob.nAngles,ii)
                imagesc(beam), colormap gray
                beamAngle = sprintf('%d^\\circ',prob.angles(ii));
                title(beamAngle,'Interpreter','tex')
                caxis([0 max(prob.x)])
                axis square
            end
            
            % Add colorbar
            pos = get(subplot(1,prob.nAngles,ii),'Position');
            cb = colorbar;
            cb.Label.String = 'Beamlet Intensity (MU)';
            set(subplot(1,prob.nAngles,ii),'Position',pos);
        end
        
        function plotBeamsPaper(prob)
            % PLOTBEAMSPAPER 
            figure()
            xRemain = prob.x;
            for ii = 1:4
                % Get beamlet intensities
                [idx,nX,nY] = FluenceMapOpt.getBeamCoords(prob.angles(ii));
                xCurrent = xRemain(1:length(idx));
                xRemain = xRemain(length(idx)+1:end);
                beam = zeros(nX,nY);
                beam(idx) = xCurrent;
                beam = beam';
                
                % Plot beamlet intensities
                subplot(2,2,ii)
                imagesc(beam), colormap gray
                set(gca,'YDir','normal','XTick',[],'YTick',[])
                caxis([0 max(prob.x)])
                axis square
            end
            
            % Positioning
            h = gcf;
            pos = h.Position;
            h.Position = [pos(1) pos(2) pos(3) pos(3)];
            a = subplot(2,2,1);
            a.Position = [0.1 0.5 0.3 0.3];
            b = subplot(2,2,2);
            b.Position = [0.45 0.5 0.3 0.3];
            c = subplot(2,2,3);
            c.Position = [0.1 0.15 0.3 0.3];
            d = subplot(2,2,4);
            d.Position = [0.45 0.15 0.3 0.3];
            
            % Colorbar
            e = colorbar('southoutside','Ticks',0:1000:3000,...
                'TickLabels',{},'LineWidth',2);
            e.Position = [0.1 0.077 0.65 0.02700];
        end
        
        function plotDose(prob,threshold)
            % PLOTDOSE Plot dose with slider for z position.
            figure()
            warning('off','MATLAB:contour:ConstantData');
            hax = axes('Units','pixels');
            
            % Plot dose at z = 50
            z = 50;
            dose = reshape(prob.D*prob.x,184,184,90);
            if nargin > 1
                idxZero = dose(:,:,:) == 0;
                dose = 2.0*(dose > threshold) - 1;
                dose(idxZero) = 0;
            end
            imagesc(dose(:,:,z),'AlphaData',dose(:,:,z)~=0), hold on
            
            % Plot body structure outlines
            for ii = 1:length(prob.mask)
                contour(prob.mask{ii}(:,:,z),1,'k');
            end
            
            % Annotations
            title(sprintf('z = %d',z))
            if nargin == 1
                caxis([min(dose(:)) max(dose(:))]);
                cb = colorbar;
                cb.Label.String = 'Dose (Gy)';
                threshold = false;
            end
            axis equal
            axis off
            hold off
            
            % Add slider
            uicontrol('Style','slider',...
                'Min',1,'Max',90,'Value',z,...
                'Position',[200 20 120 20],...
                'Callback',{@prob.updateZ,hax,dose,threshold}); 
        end
        
        function plotDosePaper(prob)
            % PLOTDOSEPAPER Plot dose at z = 50.
            figure()
            idx1 = 40:126;
            idx2 = 23:152;
            
            % Get CT slice
            ct = dicomread('Prostate_Dicom/CT.2.16.840.1.113662.2.12.0.3173.1271873797.276');
            ct = double(imresize(ct,[184,184]));
            ct50 = ct(idx1,idx2);
            ctShift = ct50 - min(ct50(:));
            ctShiftScale = ctShift/max(ctShift(:));
            CT50 = repmat(ctShiftScale,[1 1 3]);
            
            % Get Dose
            Dose = reshape(prob.D*prob.x,184,184,90);
            Dose50 = Dose(idx1,idx2,50);
            
            % Plot CT
            body50 = prob.mask{end}(idx1,idx2,50);
            imagesc(CT50), hold on
            
            % Plot dose
            imagesc(Dose50,'AlphaData',0.3*(body50~=0))
            contour(Dose50,0:10:100,'LineWidth',2);
            
            % Plot organ contours
            for i = 1:length(prob.mask)-1
               contour(prob.mask{i}(idx1,idx2,50),1,'k','LineWidth',2); 
            end
            
            % Annotations
            caxis([min(Dose50(:)) max(Dose50(:))]);
            axis equal
            axis off
            
            % Colorbar
            colorbar('southoutside','Ticks',0:20:100,'TickLabels',{},...
                'LineWidth',2)
        end
        
        function plotObjPaper(prob)
            % PLOTOBJPAPER Plot objective function values.
           
            % Objective function
            figure(1)
            subplot(3,1,1)
            plot(0:prob.nIter,prob.obj(1:prob.nIter+1),'LineWidth',2)
            FluenceMapOpt.adjustAxis(gca)
            
            % Objective terms
            for ii = 1:prob.nStructs
                for jj = 1:length(prob.structs{ii}.terms)
                    figure(1)
                    subplot(3,1,ii+1)
                    plot(0:prob.nIter,...
                        prob.structs{ii}.terms{jj}.obj(1:prob.nIter+1),...
                        'LineWidth',2);
                    FluenceMapOpt.adjustAxis(gca)
            
                    % Voxels under or over dose constraints
                    if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        figure(2), hold on
                        subplot(2,1,2)
                        plot(0:prob.nIter,...
                            prob.structs{ii}.terms{jj}.dPos(1:prob.nIter+1),...
                            'LineWidth',2);
                        FluenceMapOpt.adjustAxis(gca)
                        set(gca,'YTick',35:10:70);
                    end
                end
            end
            
            figure(2)
            subplot(2,1,1)
            plot(1:prob.nIter,prob.wDiff(1:prob.nIter),'LineWidth',2)
            FluenceMapOpt.adjustAxis(gca)
            set(gca,'YTick',0:0.05:0.15);
        end   
        
        function saveResults(prob,fileName)
            % SAVERESULTS Save current state and results.
            results = struct('structs',prob.structs,...
                'angles',prob.angles,...
                'lambda',prob.lambda,...
                'overlap',prob.overlap,...
                'x0',prob.x0,...
                'x',prob.x,...
                'obj',prob.obj,...
                'wDiff',prob.wDiff,...
                'nIter',prob.nIter,...
                'time',prob.time,...
                'tol',prob.tol,...
                'maxIter',prob.maxIter);
            save(fileName,'results');
        end
    end
    
    methods (Hidden)
        function setInputVars(prob,args)
            % SETINPUTVARS Set input variables.
            for ii = 1:length(args)/2
                prob.(args{2*ii-1}) = args{2*ii};
            end
        end 
        
        function [A,H,lb,ub] = getA(prob,type)
            % GETA Get stacked beamlet-to-voxel matrix, Hessian, and 
            %   beamlet lower and upper bounds.
            %
            %   Get output for structures and terms specified by type:
            %       * 'full': All structures and terms
            %       * 'unif': Structures and terms with uniform targets
            %       * 'slack': Columns added for slack variables
            matFull = strcmp(type,'full');
            matUnif = strcmp(type,'unif');
            matSlack = strcmp(type,'slack');

            % Add terms
            A = [];
            for ii = 1:prob.nStructs
                nVoxels = prob.structs{ii}.nVoxels;
                for jj = 1:prob.structs{ii}.nTerms
                    termUnif = strcmp(prob.structs{ii}.terms{jj}.type,'unif');
                    if (matFull || matSlack) || (matUnif && termUnif)
                        weight = prob.structs{ii}.terms{jj}.weight;
                        temp = sqrt(weight/nVoxels)*prob.structs{ii}.A;
                        if matSlack && ~termUnif
                            if strcmp(prob.structs{ii}.terms{jj}.type,'ldvc')
                                temp = -temp;
                            end
                        end
                        A = [A; temp];
                    end
                end
            end

            % Add regularization
            if prob.lambda > 0
                A = [A; sqrt(prob.lambda)*eye(prob.nBeamlets)];
            end
            
            % Add columns for slack variables
            if matSlack
                prevVoxels = 1;
                for ii = 1:prob.nStructs
                    nVoxels = prob.structs{ii}.nVoxels;
                    for jj = 1:prob.structs{ii}.nTerms
                        weight = prob.structs{ii}.terms{jj}.weight;
                        if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                            newCol = zeros(size(A,1),nVoxels);
                            eyeTerm = sqrt(weight/nVoxels)*eye(nVoxels);
                            newCol(prevVoxels:prevVoxels+nVoxels-1,:) = eyeTerm;
                            A = [A newCol];
                        end
                        prevVoxels = prevVoxels + nVoxels;
                    end
                end
            end

            % Create Hessian and beamlet bounds
            H = A'*A;
            lb = zeros(size(A,2),1);
            ub = inf(size(A,2),1);
        end
        
        function d = getd(prob,type)
            % GETD Get stacked dose vector.
            %
            %   Get output for structurs and terms specified by type:
            %       * 'full': All structures and terms
            %       * 'unif': Structures and terms with uniform targets
            %       * 'slack': Use u vectors rather than d vectors for dvcs
            vecFull = strcmp(type,'full');
            vecUnif = strcmp(type,'unif');
            vecSlack = strcmp(type,'slack');

            % Add terms
            d = [];
            for ii = 1:prob.nStructs
                nVoxels = prob.structs{ii}.nVoxels;
                for jj = 1:prob.structs{ii}.nTerms
                    termUnif = strcmp(prob.structs{ii}.terms{jj}.type,'unif');
                    if (vecFull || vecSlack) || (vecUnif && termUnif)
                        weight = prob.structs{ii}.terms{jj}.weight;
                        termD = prob.structs{ii}.terms{jj}.d;
                        temp = sqrt(weight/nVoxels)*termD;
                        if ~termUnif
                            if vecSlack
                                termU = prob.structs{ii}.terms{jj}.u;
                                if strcmp(prob.structs{ii}.terms{jj}.type,'udvc')
                                    temp = sqrt(weight/nVoxels)*termU;
                                else
                                    temp = -sqrt(weight/nVoxels)*termU;
                                end
                            else
                                termW = prob.structs{ii}.terms{jj}.w;
                                temp = temp + sqrt(weight/nVoxels)*termW;
                            end
                        end
                        d = [d; temp];
                    end
                end
            end

            % Add regularization
            if prob.lambda > 0
                d = [d; zeros(prob.nBeamlets,1)];
            end
        end
        
        function x = projX(prob,type,y0)
            % PROJX Solve non-negative least-squares problem for beamlets.
            %
            %   Solve nnls problem specified by type:
            %       * 'full': Full problem
            %       * 'unif': Initialization problem
            %       * 'slack': Problem with slack variables
            if strcmp(type,'slack')
                A = prob.As;
                H = prob.Hs;
                d = prob.getd('slack');
                x0 = y0;
                lb = prob.lbs;
                ub = prob.ubs;
            else
                if strcmp(type,'full')
                    A = prob.A;
                    H = prob.H;
                    d = prob.getd('full');
                    x0 = prob.x;
                elseif strcmp(type,'unif')
                    A = prob.Au;
                    H = prob.Hu;
                    d = prob.du;
                    x0 = zeros(prob.nBeamlets,1);
                end
                lb = prob.lb;
                ub = prob.ub;
            end
            f = -A'*d;
            if strcmp(prob.nnls,'quadprog')   
                options = optimoptions(@quadprog,'Display','off');
                x = quadprog(H,f,[],[],[],[],lb,ub,[],options);
            else
                func = @(x)FluenceMapOpt.quadObj(x,H,f);
                options.verbose = 0;
                options.method = 'newton';
                x = minConf_TMP(func,x0,lb,ub,options);
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
                        % PTV initialization
                        s = strcmp(prob.structs{ii}.terms{jj}.type,'ldvc');
                        initDose = prob.structs{ii}.A*prob.x0;
                        res = (-1)^s*(initDose - prob.structs{ii}.terms{jj}.d);
                        k = prob.structs{ii}.terms{jj}.k;
                        w = FluenceMapOpt.projW(res,k);
                        prob.structs{ii}.terms{jj}.w = w;
                        
                        % DVC paper initialization
                        % d = prob.structs{ii}.terms{jj}.d;
                        % prob.structs{ii}.terms{jj}.w = zeros(size(d));
                        
                        % OAR initialization (x0 = 0)
                        % d = prob.structs{ii}.terms{jj}.d;
                        % prob.structs{ii}.terms{jj}.w = -d;
                    end
                end
            end
        end
        
        function initU(prob)
            % INITU Initialize u vectors for dose-volume constraint terms.
            for ii = 1:prob.nStructs
                for jj = 1:prob.structs{ii}.nTerms
                    if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        % PTV initialization
                        % initDose = prob.structs{ii}.A*prob.x0;
                        % u = FluenceMapOpt.projU(0*initDose,initDose,ii,jj);
                        % prob.structs{ii}.terms{jj}.u = u;
                        
                        % DVC initialization
                        d = prob.structs{ii}.terms{jj}.d;
                        prob.structs{ii}.terms{jj}.u = d;
                        
                        % OAR initialization
                        % d = prob.structs{ii}.terms{jj}.d;
                        % prob.structs{ii}.terms{jj}.u = zeros(size(d));
                    end
                end
            end
        end
        
        function initObj(prob)
            % INITOBJ Initialize objective function values
            zeroVec = zeros(1,prob.maxIter+1);
            prob.obj = zeroVec;
            prob.wDiff = zeroVec(1:end-1);
            for ii = 1:prob.nStructs
                for jj = prob.structs{ii}.nTerms
                    prob.structs{ii}.terms{jj}.obj = zeroVec;
                    if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        prob.structs{ii}.terms{jj}.dPos = zeroVec;
                        prob.structs{ii}.terms{jj}.wPos = zeroVec;
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
                        dPos = 100*sum((-1)^s*res > 0)/nVoxels;
                        prob.structs{ii}.terms{jj}.dPos(iter+1) = dPos;
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
        
        function obj = getObj(prob,type)
            % GETOBJ Get objective value.
            if nargin == 1
                type = 'full';
            end
            if strcmp(type,'full')
                d = prob.getd('full');
                obj = 1/2*norm(prob.A*prob.x - d)^2;
            elseif strcmp(type,'unif')
                obj = 1/2*norm(prob.Au*prob.x - prob.du)^2;
            elseif strcmp(type,'dvc')
                obj = 0;
                for ii = 1:prob.nStructs
                    nVoxels = prob.structs{ii}.nVoxels;
                    for jj = 1:prob.structs{ii}.nTerms
                        if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                            weight = prob.structs{ii}.terms{jj}.weight;
                            objTerm = prob.getOarDiff(ii,jj);
                            obj = obj + weight*objTerm^2/(2*nVoxels);
                        end
                    end
                end
            else
                disp('Invalid objective type.');
                obj = -1;
            end
        end
        
        function diff = getOarDiff(prob,ii,jj,type)
            % GETOARDIFF Get OAR term ||w - (Ax - d)||_2.
            if nargin == 3
                type = 2;
            end
            Ax = prob.structs{ii}.A*prob.x;
            w = prob.structs{ii}.terms{jj}.w;
            d = prob.structs{ii}.terms{jj}.d;
            diff = norm(w - (Ax - d),type);
        end
        
        function wDiff = updateW(prob,ii,jj)
            % UPDATEW Update proxy variable.
            
            % Grab variables
            s = strcmp(prob.structs{ii}.terms{jj}.type,'ldvc');
            k = prob.structs{ii}.terms{jj}.k;
            step = prob.structs{ii}.terms{jj}.step;
            weight = prob.structs{ii}.terms{jj}.weight;
            nVoxels = prob.structs{ii}.nVoxels;
            coeff = step*weight/nVoxels;
            
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
        
        function uDiff = updateU(prob,ii,jj)
            % UPDATE U Update dose variables.
            Ax = prob.structs{ii}.A*prob.x;
            uPrev = prob.structs{ii}.terms{jj}.u;
            d = prob.structs{ii}.terms{jj}.d(1);
            n = prob.structs{ii}.nVoxels - prob.structs{ii}.terms{jj}.k;
            type = prob.structs{ii}.terms{jj}.type;
            uProj = FluenceMapOpt.projU(uPrev,Ax,d,n,type);
            step = prob.structs{ii}.terms{jj}.step;
            uDiff = norm(uPrev - uProj)/step;
            prob.structs{ii}.terms{jj}.u = uProj;
        end
        
        function [Ac,dc] = getConstraints(prob,x)
            % GETCONSTRAINTS Get stacked dose-volume constraints.
            Ac = [];
            dc = [];
            for ii = 1:prob.nStructs
                for jj = 1:prob.structs{ii}.nTerms
                    if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        [At,dt] = prob.getTermConstraint(x,ii,jj);
                        Ac = [Ac; At];
                        dc = [dc; dt];
                    end
                end
            end 
        end
        
        function [At,dt] = getTermConstraint(prob,x,ii,jj)
            % GETTERMCONSTRAINTS Get term dose-volume constraint.
            nVoxels = prob.structs{ii}.nVoxels;
            k = prob.structs{ii}.terms{jj}.k;
            s = strcmp(prob.structs{ii}.terms{jj}.type,'ldvc');
            [~,idxSort] = sort((-1)^s*prob.structs{ii}.A*x);
            At = (-1)^s*prob.structs{ii}.A(idxSort(1:nVoxels-k),:);
            dt = (-1)^s*prob.structs{ii}.terms{jj}.d(idxSort(1:nVoxels-k));
        end
        
        function grad = getIterGrad(prob,x)
            % GETITERGRAD Get gradient for iterative method.
            grad = prob.lambda*x;
            for ii = 1:prob.nStructs
                At = prob.structs{ii}.A;
                for jj = 1:prob.structs{ii}.nTerms
                    if strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        dt = prob.structs{ii}.terms{jj}.d;
                        grad = grad + At'*(At*x - dt);
                    else
                        dt = prob.getIterD(x,ii,jj);
                        if strcmp(prob.structs{ii}.terms{jj}.type,'udvc')
                            res = (At*x - dt).*(At*x > dt);
                        else
                            res = (At*x - dt).*(At*x < dt);
                        end
                        grad = grad + At'*res;
                    end
                end
            end
        end
        
        function obj = getIterObj(prob,x)
            % GETIETEROBJ Get objective for iterative method.
            obj = 1/2*norm(prob.Au*x - prob.du)^2;
            for ii = 1:prob.nStructs
                At = prob.structs{ii}.A;
                nVoxels = prob.structs{ii}.nVoxels;
                for jj = 1:prob.structs{ii}.nTerms
                    if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                        dt = prob.getIterD(x,ii,jj);
                        weight = prob.structs{ii}.terms{jj}.weight;
                        if strcmp(prob.structs{ii}.terms{jj}.type,'udvc')
                            res = (At*x - dt).*(At*x > dt);
                        else
                            res = (At*x - dt).*(At*x < dt);
                        end
                        termObj = weight/(2*nVoxels)*norm(res)^2;
                        obj = obj + termObj;
                    end
                end
            end
        end
        
        function dt = getIterD(prob,x,ii,jj)
            % GETITERD Get maximum doses for iterative method.
            n = prob.structs{ii}.nVoxels - prob.structs{ii}.terms{jj}.k;
            if strcmp(prob.structs{ii}.terms{jj}.type,'udvc')
                [~,idxSort] = sort(prob.structs{ii}.A*x,'ascend');
                dt = prob.structs{ii}.terms{jj}.d;
                dt(idxSort(n:end)) = 1e6;
            else
                [~,idxSort] = sort(prob.structs{ii}.A*x,'descend');
                dt = prob.structs{ii}.terms{jj}.d;
                dt(idxSort(n:end)) = -1e6;
            end
        end
        
        function [doses,dvh] = calcDVH(prob,x)
            % CALCDVH Calculate dose-volume histograms.
            nPoints = 1000;
            doses = linspace(0,100,nPoints);
            dvh = zeros(prob.nStructs,nPoints);
            for ii = 1:prob.nStructs
                dose = prob.structs{ii}.A*x;
                nVoxels = prob.structs{ii}.nVoxels;
                for jj = 1:nPoints
                   dvh(ii,jj) = 100*sum(dose >= doses(jj))/nVoxels;
                end
            end       
        end
        
        function plotConstraints(prob,ii,jj)
            % PLOTCONSTRAINTS Plot constraint on dose-volume histogram.
            isUnif = strcmp(prob.structs{ii}.terms{jj}.type,'unif');
            if isUnif || prob.structs{ii}.terms{jj}.percent > 0
                % Get vertical coordinates of targets/constraints
                if isUnif
                    percent = [0 100 100];
                else
                    percent = zeros(1,3);
                    constraint = prob.structs{ii}.terms{jj}.percent;
                    if strcmp(prob.structs{ii}.terms{jj}.type,'ldvc')
                        constraint = 100 - constraint;
                    end
                    percent(2:3) = constraint;
                end
                % Get horizontal coordinates of targets/constraints
                dose = zeros(1,3);
                dose(1:2) = prob.structs{ii}.terms{jj}.dose;
                plot(dose,percent,':','Color',[0.4,0.4,0.4],...
                    'LineWidth',2)
            end
        end
        
        function updateZ(prob,hObj,~,~,dose,threshold)
            % UPDATEZ Callback function for plotDose() slider.
            
            % Plot dose at current z value
            z = round(get(hObj,'Value'));
            imagesc(dose(:,:,z),'AlphaData',dose(:,:,z)~=0), hold on
            
            % Plot body structure outlines
            for ii = 1:length(prob.mask)
                contour(prob.mask{ii}(:,:,z),1,'k');
            end
            
            % Annotations
            title(sprintf('z = %d',z))
            if ~threshold
                caxis([min(dose(:)) max(dose(:))]);
                cb = colorbar;
                cb.Label.String = 'Dose (Gy)';
            end
            axis equal
            axis off
            hold off
        end
        
        function dose = getDose(prob,ii,x)
            % GETDOSE Get total dose to structure.
            if nargin == 2
                x = prob.x;
            end
            dose = sum(prob.structs{ii}.A*x);
        end
                   
        function p = getPercent(prob,ii,jj,d,x)
            % GETPERCENT Get % of voxels violating dose-volume constraint.
            if nargin == 4
                x = prob.x;
            end
            Ax = prob.structs{ii}.A*x;
            %d = prob.structs{ii}.terms{jj}.dose;
            n = prob.structs{ii}.nVoxels;
            if strcmp(prob.structs{ii}.terms{jj}.type,'udvc')
                p = 100*sum(Ax > d)/n;
            else
                p = 100*sum(Ax < d)/n;
            end
        end
        
        function doseP = getPercentile(prob,ii,p,x)
            % GETPERCENTILE Get dose at pth percentile.
            if nargin == 3
                x = prob.x;
            end
            dose = prob.structs{ii}.A*x;
            idx = floor((1-p)*length(dose));
            dose_sort = sort(dose);
            doseP = dose_sort(idx);
        end
        
        function area = getArea(prob,ii,x)
            % GETAREA Get area under dose-volume histogram curve.
            if nargin < 3
                x = prob.x;
            end
            dose = prob.structs{ii}.A*x;
            func = @(d)100*sum(dose >= d)/prob.structs{ii}.nVoxels;
            area = integral(func,0,max(dose));            
        end     
        
        function updateResults(prob,results)
           % UPDATERESULTS Update convergence results.
           %    Used for calcBeamsContinue and calcBeamsConsolidate.
            if isfield(results,'obj')
                prob.obj = results.obj;
            end
            if isfield(results,'wDiff')
                prob.wDiff = results.wDiff;
            end
            if isfield(results,'nIter')
                prob.nIter = results.nIter;
            end
            if isfield(results,'time')
                prob.time = results.time;
            end
        end
    end
    
    methods (Hidden, Static)
        function [D,nBeamlets] = getD(angles)
            % GETD Get full beamlet-to-voxel matrix.
            temp = [];
            for ii = angles
                load(['PROSTATE/Gantry' int2str(ii) '_Couch0_D.mat']);
                temp = [temp D];
            end
            D = temp;
            nBeamlets = size(D,2);
        end

        function [structs,nDVC] = getStructVars(structs,nStructs,overlap,D)
            % GETSTRUCTVARS Get structure-specific variables.
            nDVC = 0;
            vPrev = [];
            for ii = 1:nStructs
                load(['PROSTATE/' structs{ii}.name '_VOILIST.mat']);
                if ~overlap
                    [v,vPrev] = FluenceMapOpt.removeOverlap(v,vPrev); 
                end
                structs{ii}.A = D(v,:);
                structs{ii}.nVoxels = length(v);
                structs{ii}.nTerms = length(structs{ii}.terms);
                structs{ii}.terms = FluenceMapOpt.getTermVars(structs{ii}.terms,...
                    structs{ii}.nTerms,structs{ii}.nVoxels);
                for jj = 1:structs{ii}.nTerms
                   if ~strcmp(structs{ii}.terms{jj}.type,'unif')
                       nDVC = nDVC + 1;
                   end
                end
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
               load(['PROSTATE/' names{ii} '_VOILIST.mat']);
               if ~overlap
                   [v,vPrev] = FluenceMapOpt.removeOverlap(v,vPrev); 
               end
               mask{ii} = FluenceMapOpt.getMask(v);
            end
            if ~any(strcmp(names,'BODY'))
                load('PROSTATE/BODY_VOILIST.mat');
                mask{ii+1} = FluenceMapOpt.getMask(v);
            end
        end

        function mask = getMask(v)
            % GETMASK Get body structure contour for one organ.
            mask = zeros(184*184*90,1);
            mask(v) = 1;
            mask = reshape(mask,184,184,90);
        end
        
        function [fVal,gVal,hVal] = quadObj(x,H,f)
           % Objective function for non-negative least-squares problem.
           Hx = H*x;
           fVal = x'*(0.5*Hx + f);
           gVal = Hx + f;
           hVal = H;
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
        
        function u = projU(uPrev,Ax,dose,nVoxels,type)
            % PROJU Project u onto set D_v^k for slack model. 
            %
            %   udvc: D_v^k = {u : u >= uPrev & #(u <= dose) >= nVoxel}
            %   ldvc: D_v^k = {u : u <= uPrev & #(u >= dose) >= nVoxel}
            count = 0;
            if strcmp(type,'udvc')
                u = max(uPrev,Ax);
                [~,idx] = sort(u,'ascend');
                for ii = 1:length(u)
                    if u(idx(ii)) <= dose
                        count = count + 1;
                    else
                        if uPrev(idx(ii)) <= dose
                            u(idx(ii)) = dose;
                            count = count + 1;
                        end
                    end
                    if count >= nVoxels
                        break
                    end
                end
            else
                u = min(uPrev,Ax);
                [~,idx] = sort(u,'descend');
                for ii = 1:length(u)
                    if u(idx(ii)) >= dose
                        count = count + 1;
                    else
                        if uPrev(idx(ii)) >= dose
                            u(idx(ii)) = dose;
                            count = count + 1;
                        end
                    end
                    if count >= nVoxels
                        break
                    end
                end
            end
        end
        
        function [idx,nX,nY] = getBeamCoords(angle)
            % GETBEAMCOORDS Get beamlet coordinates.
            load(['PROSTATE/Gantry' int2str(angle) '_Couch0_BEAMINFO.mat']);
            xIdx = x - min(x) + 1;
            yIdx = y - min(y) + 1;
            nX = max(xIdx);
            nY = max(yIdx);
            idx = sub2ind([nX,nY],xIdx,yIdx);
        end
        
        function adjustAxis(g)
            % ADJUSTAXIS Readjust axes limits for objective plot.
            axis tight
            yVals = g.YLim;
            yPad = 0.1*(yVals(2) - yVals(1));
            g.YLim = [yVals(1)-yPad yVals(2)+yPad];
            g.XLim = [0 42];
            g.XTick = 0:6:43;
            g.XTickLabels = {};
            g.YTickLabels = {};
            g.LineWidth = 2;      
        end
    end
end
