function prob = calcBeamsConsolidate(prob,print)
    % CALCBEAMSCONSOLIDATE Consolidate all dose-volume terms for each
    % structure into one objective term.
    
    % Initialization
    tic;
    prob.x = prob.x0;
    [A,H,nDVC] = getA(prob);
    [yMat,~] = projY(prob);
    obj = getObj(prob,A,yMat);
    wDiff = [];
    if print
        fprintf('iter: %d, obj: %7.4e\n',0,obj);
    end
    
    % Fluence map optimization
    for kk = 1:prob.maxIter
        prob.x = projX(prob,A,H,yMat);
        [yMat,yDiff] = projY(prob,yMat,nDVC);
        wDiff = [wDiff yDiff];
        obj = [obj getObj(prob,A,yMat)];
        if print
            fprintf('iter: %d, obj: %7.4e, yDiff: %7.4e\n',kk,...
                obj(end),yDiff);
        end
        if yDiff <= prob.tol
            break
        end
    end
    prob.x = projX(prob,A,H,yMat);
end

function [A,H,nDVC] = getA(prob)
    % Assumes at most one unif term per structure
    A = [];
    nDVC = 0;
    for ii = 1:prob.nStructs
        structA = prob.structs{ii}.A;
        structVoxels = prob.structs{ii}.nVoxels;
        addDvcTerm = true;
        for jj = 1:prob.structs{ii}.nTerms
            termWeight = prob.structs{ii}.terms{jj}.weight;
            termA = sqrt(termWeight/structVoxels)*structA;
            if strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                A = [A; termA];
            elseif addDvcTerm
                A = [A; termA];
                nDVC = nDVC + 1;
                addDvcTerm = false;
            end 
        end
    end
    A = [A; sqrt(prob.lambda)*eye(prob.nBeamlets)];
    H = A'*A;
end

function [yMat,yDiff] = projY(prob,yMat,nDVC)
    % Assumes at least one DVC term per structure
    yDiff = 0;
    for ii = 1:prob.nStructs
        y = prob.structs{ii}.A*prob.x;
        for jj = 1:prob.structs{ii}.nTerms
            [~,idxAscend] = sort(y,'ascend');
            [~,idxDescend] = sort(y,'descend');
            d = prob.structs{ii}.terms{jj}.d(1);
            if strcmp(prob.structs{ii}.terms{jj}.type,'udvc')
                k = prob.structs{ii}.terms{jj}.k;
                y(idxDescend(k+1:end)) = min(y(idxDescend(k+1:end)),d);
            elseif strcmp(prob.structs{ii}.terms{jj}.type,'ldvc')
                k = prob.structs{ii}.terms{jj}.k;
                y(idxAscend(k+1:end)) = max(y(idxAscend(k+1:end)),d);
            end
        end
        if nargin > 1
            step = prob.structs{ii}.terms{jj}.step;
            yDiff = yDiff + norm(yMat{ii} - y)/(step*nDVC);
        end
        yMat{ii} = y;
    end
end

function x = projX(prob,A,H,yMat)
    d = getd(prob,yMat);
    f = -A'*d;
    func = @(x)FluenceMapOpt.quadObj(x,H,f);
    options.verbose = 0;
    options.method = 'newton';
    x = minConf_TMP(func,prob.x,zeros(prob.nBeamlets,1),...
        inf(prob.nBeamlets,1),options);
end

function d = getd(prob,yMat)
    d = [];
    for ii = 1:prob.nStructs
        structVoxels = prob.structs{ii}.nVoxels;
        addDvcTerm = true;
        for jj = 1:prob.structs{ii}.nTerms
            termWeight = prob.structs{ii}.terms{jj}.weight;
            termD = prob.structs{ii}.terms{jj}.d;
            if strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                termD = sqrt(termWeight/structVoxels)*termD;
                d = [d; termD];
            elseif addDvcTerm
                termD = sqrt(termWeight/structVoxels)*yMat{ii};
                d = [d; termD];
                addDvcTerm = false;
            end
        end
    end
    d = [d; zeros(prob.nBeamlets,1)];
end

function obj = getObj(prob,A,yMat)
    d = getd(prob,yMat);
    obj = 1/2*norm(A*prob.x-d)^2;
end

function diff = getDiff(prob,yMat)
    diff = 0;
    for ii = 1:prob.nStructs
        structA = prob.structs{ii}.A;
        addDvcTerm = true;
        for jj = 1:prob.structs{ii}.nTerms
            if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif') && addDvcTerm
                termDiff = norm(structA*prob.x - yMat{ii},inf);
                diff = max(diff,termDiff);
                addDvcTerm = false;
            end
        end
    end
end
