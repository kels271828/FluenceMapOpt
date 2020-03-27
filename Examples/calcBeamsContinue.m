function prob = calcBeamsContinue(prob,gamma,sigma,tol,maxIter,...
    print,consolidate)
    % CALCBEAMSCONTINUE Approach inspired by Lu paper.
    
    % Initialization
    t1 = clock;
    prob.initProb(false);
    obj = prob.getObj();
    oDiff = getDiff(prob);
    if print
        fprintf('iter: 0, obj: %7.4e, oDiff: %7.4e\n',obj,oDiff);
    end
    
    % Continuation
    for kk = 1:maxIter
        % Solve approximate solution
        if consolidate
           [prob,objTemp,oDiffTemp] = calcBeamsConsolidate(prob,false);
           obj = [obj objTemp];
           oDiff = [oDiff oDiffTemp];
        else
            prob.calcBeams(false);
            obj = [obj prob.getObj()];
            oDiff = [oDiff getDiff(prob)];
        end
            
        % Print status
        if print
            fprintf('iter: %d, obj: %7.4e, oDiff: %7.4e\n',kk,obj(end),oDiff(end));
        end
        
        % Check convergence
        if oDiff(end) <= tol
            break
        end
        
        % Update parameters
        for ii = 1:prob.nStructs
            nVoxels = prob.structs{ii}.nVoxels;
            for jj = 1:prob.structs{ii}.nTerms
                if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                    oldWeight = prob.structs{ii}.terms{jj}.weight;
                    newWeight = sigma*oldWeight;
                    prob.structs{ii}.terms{jj}.weight = sigma*oldWeight;
                    prob.structs{ii}.terms{jj}.step = nVoxels/newWeight;
                end
            end
        end
        [prob.A,prob.H,~,~] = prob.getA('full');
        prob.tol = gamma*prob.tol;
        prob.x0 = prob.x;  
    end
    prob.nIter = kk;
    prob.obj = obj;
    prob.wDiff = oDiff;
    t2 = clock;
    prob.time = etime(t2,t1);
end

function diff = getDiff(prob)
    diff = 0;
    for ii = 1:prob.nStructs
        for jj = 1:prob.structs{ii}.nTerms
            if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                diff = max(diff,prob.getOarDiff(ii,jj,Inf));
            end
        end
    end
end