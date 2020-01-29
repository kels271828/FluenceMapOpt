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
        names    % Body structure names
        mask     % Body structure contours
        D        % Full beamlet-to-voxel matrix
        
        nStructs % Number of structures
        nAngles  % Number of angles
        nBeamlts % Bumber of beamlets
        
        
        
        A        % Stacked beamlet-to-voxel matrix
        d        % Stacked dose vector
        H        % Hessian of objective function
        A_unif   % Stacked beamlet-to-voxel matrix for uniform dose targets
        d_unif   % Stacked dose vector for uniform dose targets
        H_unif   % Stacked Hessian matrix for uniform dose targets
        lb       % Lower bound for beamlets
        ub       % Upper bound for beamlets
        
        % Solution variables
        x        % Final beamlet intensities
        obj      % Objective function values
        err      % Error values between w vectors
        nIter    % Number of iterations used
    end
    
    methods
        function prob = FluenceMapOpt(structs,varargin)
        % FLUENCEMAPOPT Initialize problem.
        
            % Set input variables
            if nargin >= 2, prob.angles = varargin{2}; end
            if nargin >= 3, prob.lambda = varargin{3}; end
            if nargin >= 4, prob.maxIter = varargin{4}; end
            if nargin >= 5, prob.tol = varargin{5}; end
            if nargin == 6, prob.overlap = varargin{6}; end
            
            % Compute internal variables
            prob.D = getD(prob.angles);
            prob.structs = getStructVars(structs,prob.overlap,prob.D);
            prob.names = getNames(prob.structs);
            prob.mask = getMaskStruct(prob.names,prob.overlap);
               
            prob.nStructs = length(prob.structs);
            prob.nAngles = length(prob.angles);
            prob.nBeamlts = size(prob.D,2);
            
            prob.getA('full'); 
            
            % Compute initial beamlets
            if nargin >= 1
                prob.xInit = varargin{1};
            else
                prob.xInit = prob.getInitX();
            end
            
            % Set solution variables if given
            if flag && isfield(pars,'x')
                f.x = pars.x;
            else
                f.x = f.xInit;
            end
            if flag && isfield(pars,'obj')
                f.obj = pars.obj;
            end
            if flag && isfield(pars,'nIter')
                f.nIter = pars.nIter;
            end
        end
        
        % Get stacked A matrix.
        function getA(f,type)
            
            A = [];
            for i = 1:f.nStructs
                for j = 1:f.structs{i}.nTerms
                    temp = sqrt(f.structs{i}.terms{j}.weight/f.structs{i}.nVoxels)*f.structs{i}.A;
                    if strcmp(type,'unif')
                        if strcmp(f.structs{i}.terms{j}.type,'unif')
                            A = [A; temp];
                        end
                    else
                        A = [A; temp];
                    end
                end
            end
            if f.lambda > 0
                A = [A; sqrt(f.lambda)*eye(f.nBeamlts)];
            end
            if strcmp(type,'unif')
                f.A_unif = A;
                f.H_unif = A'*A;
            else
                f.A = A;
                f.H = A'*A;
            end
            f.lb = zeros(f.nBeamlts,1);
            f.ub = inf(f.nBeamlts,1);
        end
        
        % Get stacked d vector.
        function getd(f,type)
            
            d = [];
            for i = 1:f.nStructs
                for j = 1:f.structs{i}.nTerms
                    temp = sqrt(f.structs{i}.terms{j}.weight/f.structs{i}.nVoxels)*f.structs{i}.terms{j}.d;
                    if strcmp(type,'unif')
                        if strcmp(f.structs{i}.terms{j}.type,'unif')
                            d = [d; temp];
                        end
                    else
                        if ~strcmp(f.structs{i}.terms{j}.type,'unif')
                            temp = temp + sqrt(f.structs{i}.terms{j}.weight/f.structs{i}.nVoxels)*f.structs{i}.terms{j}.w;
                        end
                        d = [d; temp];
                    end
                end
            end
            if f.lambda > 0
                d = [d; zeros(f.nBeamlts,1)];
            end
            if strcmp(type,'unif')
                f.d_unif = d;
            else
                f.d = d;
            end
        end
       
        % Initialize x.
        function initX(f)
            
            f.getA('unif');
            f.getd('unif');

            % Solve non-negative least-squares problem for x           
            F = -f.A_unif'*f.d_unif;
            fun = @(x)f.quadObj(x,f.H_unif,F);
            options.verbose = 0;
            options.method = 'newton';
            f.xInit = minConf_TMP(fun,zeros(f.nBeamlts,1),f.lb,f.ub,options);
        end

        % Constraint generation method.
        function constGen(f)
            
            
            [A_dvc,d_dvc] = f.calcConstMatBig(f.xInit);
            f.x = lsqlin(f.A_unif,f.d_unif,A_dvc,d_dvc,[],[],f.lb,f.ub);
        end
        
        % Calculate stacked constraint matrix for dose-volume constraints.
        function [A_dvc,d_dvc] = calcConstMatBig(f,x)
           
            A_dvc = [];
            d_dvc = [];
            for i = 1:f.nStructs
                for j = 1:f.structs{i}.nTerms
                    if ~strcmp(f.structs{i}.terms{j}.type,'unif')
                        [A_temp,d_temp] = calcConstMatSmall(f,x,i,j);
                        A_dvc = [A_dvc; A_temp];
                        d_dvc = [d_dvc; d_temp];
                    end
                end
            end
        end
        
        % Calculate constraint matrix for dose-volume constraint.
        function [A_dvc,d_dvc] = calcConstMatSmall(f,x,struct,term)
            
            [~,idx] = sort(f.structs{struct}.A*x);
            if strcmp(f.structs{struct}.terms{term}.type,'udvc')
                idx_lower = f.structs{struct}.nVoxels - f.structs{struct}.terms{term}.k;
                A_dvc = f.structs{struct}.A(idx(1:idx_lower),:);
                d_dvc = f.structs{struct}.terms{term}.d(idx(1:idx_lower));
            elseif strcmp(f.structs{struct}.terms{term}.type,'ldvc')
                idx_upper = f.structs{struct}.terms{term}.k + 1;
                A_dvc = -f.structs{struct}.A(idx(idx_upper:end),:);
                d_dvc = -f.structs{struct}.terms{term}.d(idx(idx_upper:end));
            else
                A_dvc = NaN;
                d_dvc = NaN;
            end 
        end
        
        % conrad method.
        function conrad(f,slope)
            
            if ~exist('slope','var')
                slope = 1;
            end

            % Non-negative least-squares with relaxed OAR constraint
            cvx_begin quiet
                variable x1(f.nBeamlts)
                minimize( sum_square(f.A_unif*x1 - f.d_unif) )
                subject to
                    f.lb <= x1;
                    for i = 1:f.nStructs
                        for j = 1:f.structs{i}.nTerms
                            if ~strcmp(f.structs{i}.terms{j}.type,'unif')
                                A = f.structs{i}.A;
                                d = f.structs{i}.terms{j}.d;
                                if strcmp(f.structs{i}.terms{j}.type,'udvc')
                                    k = f.structs{i}.nVoxels - f.structs{i}.terms{j}.k;
                                else
                                    k = f.structs{i}.terms{j}.k;
                                end
                                sum(pos(1 + slope*(A*x1 - d))) <= k;
                            end
                        end
                    end
            cvx_end

            % Non-negative least-squares with hard OAR constraint
            cvx_begin quiet
                variable x2(f.nBeamlts)
                minimize( sum_square(f.A_unif*x2 - f.d_unif) )
                subject to
                    f.lb <= x2;
                    for i = 1:f.nStructs
                        for j = 1:f.structs{i}.nTerms
                            if ~strcmp(f.structs{i}.terms{j}.type,'unif')
                                [A,d] = calcConstMatSmall(f,x1,i,j);
                                A*x2 <= d;
                            end
                        end
                    end
            cvx_end
            
            f.x = x2;
        end
        
        % Calculate beamlet intensities.
        function calcDose(f,print)
            
            if ~exist('print','var')
                print = false;
            end
            
            % Initialize x, w, and objective values
            f.initProb(print);
            
            % Fluence map optimization
            for t = 1:f.maxIter

                % Update x
                f.projX();
                
                % Update w
                errorSum = 0;
                for i = 1:f.nStructs
                    for j = 1:f.structs{i}.nTerms
                        if ~strcmp(f.structs{i}.terms{j}.type,'unif')
                            
                            % Update wPrev
                            wPrev = f.structs{i}.terms{j}.w;
                            
                            % Update w
                            Axmd = f.structs{i}.A*f.x - f.structs{i}.terms{j}.d;
                            coeff = f.structs{i}.terms{j}.step*f.structs{i}.terms{j}.weight/f.structs{i}.nVoxels;
                            temp = f.structs{i}.terms{j}.w + coeff*(Axmd - f.structs{i}.terms{j}.w);
                            s = strcmp(f.structs{i}.terms{j}.type,'ldvc');
                            f.structs{i}.terms{j}.w = (-1)^s*f.projW((-1)^s*temp,f.structs{i}.terms{j}.k);
                            
                            % Error sum
                            errorSum = errorSum + norm(f.structs{i}.terms{j}.w - wPrev)/f.structs{i}.terms{j}.step;  
                        end
                    end
                end
                f.err(t) = errorSum;
                f.nIter = t;
                f.calcObj(t,print);

                % Stopping criteria
                if errorSum <= f.tol
                    break
                end
                
            end
        end

        % Initialize x, w, and objective values
        function initProb(f,print)

            % Initialize x and w
            f.x = f.xInit;
            for i = 1:f.nStructs
                for j = 1:f.structs{i}.nTerms
                    if ~strcmp(f.structs{i}.terms{j}.type,'unif')
                        temp = f.structs{i}.A*f.x - f.structs{i}.terms{j}.d;
                        s = strcmp(f.structs{i}.terms{j}.type,'ldvc');
                        f.structs{i}.terms{j}.w = (-1)^s*f.projW((-1)^s*temp,f.structs{i}.terms{j}.k);
                    end
                end
            end
            
            % Initialize objective function values
            f.obj = zeros(1,f.maxIter+1);
            f.err = zeros(1,f.maxIter);
            for i = 1:f.nStructs
                for j = 1:f.structs{i}.nTerms
                    f.structs{i}.terms{j}.obj = zeros(1,f.maxIter+1);
                    if ~strcmp(f.structs{i}.terms{j}.type,'unif')
                        f.structs{i}.terms{j}.vdiff = zeros(1,f.maxIter+1);
                        f.structs{i}.terms{j}.wdiff = zeros(1,f.maxIter+1);
                    end
                end
            end
            
            % Calculate and print initial objective value
            f.calcObj(0,print)
        end
        
        % Solve non-negative least-squares problem for x.
        function projX(f)
            
            f.getd('full');
            F = -f.A'*f.d;
            fun = @(x)f.quadObj(x,f.H,F);
            options.verbose = 0;
            options.method = 'newton';
            f.x = minConf_TMP(fun,f.x,f.lb,f.ub,options);
        end
        
        % Objective function for non-negative least squares problem for x.
        function [fval,gval,Hval] = quadObj(~,x,H,F)
            
            Hx = H*x;
            fval = x'*(0.5*Hx + F);
            gval = Hx + F;
            Hval = H;
        end
        
        % Keep k largest positive entries of w and set the rest to zero.
        function w = projW(~,w,k)
            
            idxPos = w > 0;
            wPos = w(idxPos);
            if sum(idxPos) > k
                [~,idxSort] = sort(wPos,'descend');
                wPos(idxSort(k+1:end)) = 0;
                w(idxPos) = wPos;
            end
        end
        
        % Calculate and print objective function value.
        function calcObj(f,iter,print)
            
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
        
        % Plot objective function values.
        function plotObj(f)
            
            figure()
            myLines = lines;
           
            % Objective function
            subplot(2,2,1)
            plot(0:f.nIter,f.obj(1:f.nIter+1),'Color',[0.5,0.5,0.5])
            xlabel('Iteration (k)')
            ylabel('Objective Value')
            
            % Convergence of w variables
            subplot(2,2,3)
            plot(1:f.nIter,f.err(1:f.nIter),'Color',[0.5,0.5,0.5])
            xlabel('Iteration (k)')
            ylabel('Convergence Criteria')
            
            for i = 1:f.nStructs
                for j = 1:length(f.structs{i}.terms)
                    % Objective function terms
                    subplot(2,2,2), hold on
                    plot(0:f.nIter,f.structs{i}.terms{j}.obj(1:f.nIter+1),'Color',myLines(i,:));
            
                    % Voxels under or over dose constraints
                    if ~strcmp(f.structs{i}.terms{j}.type,'unif')
                        subplot(2,2,4), hold on
                        plot(0:f.nIter,f.structs{i}.terms{j}.vdiff(1:f.nIter+1),'Color',myLines(i,:));
                        plot(0:f.nIter,f.structs{i}.terms{j}.wdiff(1:f.nIter+1),'--','Color',myLines(i,:));
                    end
                end
            end
            subplot(2,2,2)
            xlabel('Iteration (k)')
            ylabel('Objective Terms')
            legend(f.names,'Location','NorthEast')
            
            subplot(2,2,4)
            xlabel('Iteration (k)')
            ylabel('% Voxels Exceeding Dose')
        end
        
        % Plot objective function values (fig 8).
        function plotObjPaper(f)
            
            myLines = lines;
           
            % Objective function
            figure(1)
            subplot(3,1,1)
            plot(0:f.nIter,f.obj(1:f.nIter+1),'Color',[0.5,0.5,0.5],'LineWidth',3)
            f.adjustAxis(gca)
            
            % Objective terms
            for i = 1:f.nStructs
                for j = 1:length(f.structs{i}.terms)
                    figure(1)
                    subplot(3,1,i+1)
                    plot(0:f.nIter,f.structs{i}.terms{j}.obj(1:f.nIter+1),'Color',myLines(i,:),'LineWidth',3);
                    f.adjustAxis(gca)
            
                    % Voxels under or over dose constraints
                    if ~strcmp(f.structs{i}.terms{j}.type,'unif')
                        figure(2), hold on
                        subplot(2,1,2)
                        plot(0:f.nIter,f.structs{i}.terms{j}.vdiff(1:f.nIter+1),'Color',myLines(i,:),'LineWidth',3);
                        f.adjustAxis(gca)
                        set(gca,'YTick',52:2:56);
                    end
                end
            end
            
            figure(2)
            subplot(2,1,1)
            plot(1:f.nIter,f.err(1:f.nIter),'Color',[0.5,0.5,0.5],'LineWidth',3)
            f.adjustAxis(gca);
        end   
        
        % Readjust axes limits.
        function adjustAxis(~,g)
            
            axis tight
            yVals = g.YLim;
            yPad = 0.1*(yVals(2) - yVals(1));
            g.YLim = [yVals(1)-yPad yVals(2)+yPad];
            g.XTick = 0:50:200;
            g.XTickLabels = {};
            g.YTickLabels = {};
            g.LineWidth = 2;      
        end
        
        % Calculate and plot dose-volume histogram of solution.
        function plotDVH(f)
            
            myLines = lines;
            
            % Calculate dose-volume histograms
            doses = linspace(0,100,1000);
            dvhInit = zeros(f.nStructs,length(doses));
            dvhFinal = zeros(f.nStructs,length(doses));
            for i = 1:f.nStructs
                doseInit = f.structs{i}.A*f.xInit;
                doseFinal = f.structs{i}.A*f.x;
                for j = 1:length(doses)
                    dvhInit(i,j) = 100*sum(doseInit > doses(j))/f.structs{i}.nVoxels;
                    dvhFinal(i,j) = 100*sum(doseFinal > doses(j))/f.structs{i}.nVoxels;
                end
            end
            
            % Plot dose-volume histograms
            figure(), hold on
            
            legendHandles = [];
            legendNames = {};
            for i = 1:f.nStructs
                for j = 1:length(f.structs{i}.terms)
                    if ~strcmp(f.structs{i}.terms{j}.type,'unif') && f.structs{i}.terms{j}.percent == 0
                        plot(f.structs{i}.terms{j}.dose,0,'p','MarkerFaceColor',[0.9290 0.6940 0.1250],...
                            'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerSize',10);
                    else
                        if strcmp(f.structs{i}.terms{j}.type,'unif')
                            percent = [0 100 100];
                        elseif f.structs{i}.terms{j}.percent > 0
                            percent = zeros(1,3);
                            percent(2:3) = f.structs{i}.terms{j}.percent;
                        end
                        dose = zeros(1,3);
                        dose(1:2) = f.structs{i}.terms{j}.dose;
                        plot(dose,percent,':','Color',[0.4 0.4 0.4])
                        plot(doses,dvhInit(i,:),'--','Color',myLines(i,:))
                        if j == 1
                            lineHandle = plot(doses,dvhFinal(i,:),'Color',myLines(i,:));
                            lineName = f.structs{i}.name;
                            legendHandles = [legendHandles lineHandle];
                            legendNames = [legendNames, lineName];
                        else
                            plot(doses,dvhFinal(i,:),'Color',myLines(i,:))
                        end
                    end
                end
                
                % Annotations
                legend(legendHandles,legendNames)
                xlabel('Dose (Gy)')
                ylabel('Relative Volume (%)')
                ax = gca;
                ax.XLim = [0 doses(end)];
                ax.YLim = [0 100];
                box on
                axis square
            end
        end
        
        % Calculate and plot dose-volume histogram of solution (fig 9,11,12,13).
        function plotDVHPaper(f)
            
            myLines = lines;
            
            % Calculate dose-volume histograms
            doses = linspace(0,100,1000);
            dvhInit = zeros(f.nStructs,length(doses));
            dvhFinal = zeros(f.nStructs,length(doses));
            for i = 1:f.nStructs
                doseInit = f.structs{i}.A*f.xInit;
                doseFinal = f.structs{i}.A*f.x;
                for j = 1:length(doses)
                    dvhInit(i,j) = 100*sum(doseInit > doses(j))/f.structs{i}.nVoxels;
                    dvhFinal(i,j) = 100*sum(doseFinal > doses(j))/f.structs{i}.nVoxels;
                end
            end
            
            % Plot dose-volume histograms
            for i = 1:f.nStructs
                figure(), hold on
                for j = 1:length(f.structs{i}.terms)
                    if ~strcmp(f.structs{i}.terms{j}.type,'unif') && f.structs{i}.terms{j}.percent == 0
                        plot(f.structs{i}.terms{j}.dose,0,'p','MarkerFaceColor',[0.9290 0.6940 0.1250],...
                            'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerSize',10);
                    else
                        if strcmp(f.structs{i}.terms{j}.type,'unif')
                            percent = [0 100 100];
                        elseif f.structs{i}.terms{j}.percent > 0
                            percent = zeros(1,3);
                            percent(2:3) = f.structs{i}.terms{j}.percent;
                        end
                        dose = zeros(1,3);
                        dose(1:2) = f.structs{i}.terms{j}.dose;
                        plot(dose,percent,':','Color',[0.4 0.4 0.4],'LineWidth',3)
                        plot(doses,dvhInit(i,:),'--','LineWidth',3,'Color',myLines(i,:))
                        plot(doses,dvhFinal(i,:),'LineWidth',3,'Color',myLines(i,:))
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
        
        % Plot beamlet intensities.
        function plotBeamlets(f)
            
            figure()
            x = f.x;
            
            for i = 1:f.nAngles
                % Get x and y positions
                [linIdx,nx,ny] = f.getBeamlets(f.angles(i));
                C = zeros(nx,ny);
                C(linIdx) = 1;
                
                % Get beamlet intensities
                xTemp = x(1:length(linIdx));
                x = x(length(linIdx)+1:end);
                B = zeros(nx,ny);
                B(linIdx) = xTemp;
                B = B';
                
                % Plot beamlet intensities
                subplot(1,f.nAngles,i)
                imagesc(B), colormap gray
                title(sprintf('%d^\\circ',f.angles(i)),'Interpreter','tex')
                set(gca,'YDir','normal','XTick',[],'YTick',[])
                caxis([0 max(f.x)])
                axis square
            end
            p = get(subplot(1,f.nAngles,i),'Position');
            cb = colorbar;
            cb.Label.String = 'Beamlet Intensity (MU)';
            set(subplot(1,f.nAngles,i),'Position',p);
        end
        
        % Plot beamlet intensities (fig 10,14).
        function plotBeamletsPaper(f)
            
            figure()
            x = f.x;
            
            for i = 1:4
                % Get x and y positions
                [linIdx,nx,ny] = f.getBeamlets(f.angles(i));
                
                % Get beamlet intensities
                xTemp = x(1:length(linIdx));
                x = x(length(linIdx)+1:end);
                B = zeros(nx,ny);
                B(linIdx) = xTemp;
                B = B';
                
                % Plot beamlet intensities
                subplot(2,2,i)
                imagesc(B), colormap gray
                set(gca,'YDir','normal','XTick',[],'YTick',[])
                caxis([0 max(f.x)])
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
            e = colorbar('southoutside','Ticks',0:1000:3000,'TickLabels',{},'LineWidth',2);
            e.Position = [0.1    0.077    0.65    0.02700];
        end
        
        % Get x and y positions for beamlets.
        function [linIdx,nx,ny] = getBeamlets(~,angle)
            
            load(['Gantry' int2str(angle) '_Couch0_BEAMINFO.mat']);
            xIdx = x - min(x) + 1;
            yIdx = y - min(y) + 1;
            nx = max(xIdx);
            ny = max(yIdx);
            linIdx = sub2ind([nx,ny],xIdx,yIdx);
        end
        
        % Plot dose with slider for z position.
        function plotDose(f)
            
            figure()
            
            % Get dose
            Dose = reshape(f.D*f.x,184,184,90);
            minDose = min(Dose(:));
            maxDose = max(Dose(:));
            str = 'z = %d';
            c = [minDose maxDose];
            
            % Plot dose and contours
            warning('off','MATLAB:contour:ConstantData');
            hax = axes('Units','pixels');
            imagesc(Dose(:,:,50),'AlphaData',Dose(:,:,50)~=0), hold on
            for i = 1:length(f.mask)
                contour(f.mask{i}(:,:,50),1,'k');
            end
            
            % Annotations
            title(sprintf(str,50))
            caxis(c);
            cb = colorbar;
            cb.Label.String = 'Dose (Gy)';
            axis equal
            axis off
            hold off
            
            % Add slider
            uicontrol('Style','slider',...
                'Min',1,'Max',90,'Value',50,...
                'Position',[200 20 120 20],...
                'Callback',{@f.whichSlice,hax,Dose,str,c});
        end
        
        % Plot dose at slice 50 (fig 1,5,10,14).
        function plotDosePaper(f)
            
            figure()
            idx1 = 40:126;
            idx2 = 23:152;
            
            % Get CT slice
            ct = dicomread('CT.2.16.840.1.113662.2.12.0.3173.1271873797.276');
            ct = double(imresize(ct,[184,184]));
            ct50 = ct(idx1,idx2);
            ctShift = ct50 - min(ct50(:));
            ctShiftScale = ctShift/max(ctShift(:));
            CT50 = repmat(ctShiftScale,[1 1 3]);
            
            % Get Dose
            Dose = reshape(f.D*f.x,184,184,90);
            Dose50 = Dose(idx1,idx2,50);
            
            % Plot CT
            body50 = f.mask{end}(idx1,idx2,50);
            imagesc(CT50), hold on
            
            % Plot dose
            imagesc(Dose50,'AlphaData',0.3*(body50~=0))
            contour(Dose50,0:10:100,'LineWidth',2);
            
            % Plot organ contours
            for i = 1:length(f.mask)-1
               contour(f.mask{i}(idx1,idx2,50),1,'k','LineWidth',2); 
            end
            
            % Annotations
            caxis([min(Dose50(:)) max(Dose50(:))]);
            axis equal
            axis off
            
            % Colorbar
            colorbar('southoutside','Ticks',0:20:100,'TickLabels',{},'LineWidth',2)
        end
        
        % Callback function for plotDose() slider.
        function whichSlice(f,hObj,~,~,Dose,str,c)
            
            % Plot dose and contours
            z = round(get(hObj,'Value'));
            imagesc(Dose(:,:,z),'AlphaData',Dose(:,:,z)~=0), hold on
            for i = 1:length(f.mask)
                contour(f.mask{i}(:,:,z),1,'k');
            end
            
            % Annotations
            title(sprintf(str,z))
            caxis(c);
            cb = colorbar;
            cb.Label.String = 'Dose (Gy)';
            axis equal
            axis off   
            hold off
        end
        
        % Save results
        function saveResults(f,str)
            
            % Input parameters
            pars.structs = f.structs;
            pars.angles = f.angles;
            pars.lambda = f.lambda;
            pars.maxIter = f.maxIter;
            pars.xInit = f.xInit;
            pars.overlap = f.overlap;
            pars.tol = f.tol;
            
            % Solution variables
            pars.x = f.x;
            pars.obj = f.obj;
            pars.nIter = f.nIter;
            save(str,'pars');
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

function [v,vPrev] = removeOverlap(v,vPrev)
% REMOVEOVERLAP Remove overlapping voxels from body structure.

   if isempty(vPrev)
       vPrev = v;
   else
       v = setdiff(v,vPrev);
       vPrev = union(v,vPrev);
   end 
end

function structs = getStructVars(structs,overlap,D)
% GETSTRUCTVARS Get structure-specific variables.

    vPrev = [];
    for ii = 1:length(structs)
        load(structs{ii}.name '_VOILIST.mat']);
        if ~overlap, [v,vPrev] = removeOverlap(v,vPrev); end
        structs{ii}.A = D(v,:);
        structs{ii}.terms = getTermVars(structs{ii}.terms,length(v));
    end
end

function terms = getTermVars(terms,nVoxels)
% GETTERMVARS Get term-specific variables.

    for ii = 1:length(terms)
       terms{ii}.d = terms{ii}.dose*ones(nVoxels,1);
       terms{ii}.step = nVoxels/terms{ii}.weight;
       if ~strcmp(terms{ii}.type,'udvc')
           terms{ii}.k = floor(terms{ii}.percent*nVoxels/100);
       end
    end
end

function names = getNames(structs)
% GETNAMES Get body structure names.

    nStructs = length(structs);
    names = cell(1,nStructs);
    for ii = 1:nStructs
        names{ii} = prob.structs{ii}.name;
    end
end

function mask = getMask(v)
% GETMASK Get body structure contour for one organ.

    mask = zeros(184*184*90,1);
    mask(v) = 1;
    mask = reshape(mask,184,184,90);
end

function mask = getMaskStruct(names,overlap)
% GETMASKSTRUCT Get body structure contours for all organs.

    vPrev = [];
    for ii = 1:length(names)
       load([names{ii} '_VOILIST.mat']);
       if ~overlap, [v,vPrev] = removeOverlap(v,vPrev); end
       mask{ii} = getMask(v);
    end
    if ~any(strcmp(names,'BODY'))
        load('BODY_VOILIST.mat');
        mask{ii+1} = getMask(v);
    end
end
