function prob = calcBeamsCont(prob,structs,print)
    % CALCBEAMSCONT Approach inspired by Lu paper.
    
    % Initialization
    t1 = clock;
    prob.initProb(false);
    obj = prob.getObj();
    oDiff = getDiff(prob);
    if print
        fprintf('iter: 0, obj: %7.4e, oDiff: %7.4e\n',obj,oDiff);
    end
    
    % Continuation
    gamma = 0.5;
    sigma = 1.5;
    tol = prob.tol;
    for kk = 1:prob.maxIter
        % Solve approximate solution
        prob.calcBeams(false);
        obj = [obj prob.getObj()];
        oDiff = [oDiff getDiff(prob)];
        
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
            for jj = 1:prob.structs{ii}.nTerms
                if ~strcmp(structs{ii}.terms{jj}.type,'unif')
                    oldWeight = structs{ii}.terms{jj}.weight;
                    structs{ii}.terms{jj}.weight = sigma*oldWeight;
                end
            end
        end
        prob.tol = gamma*prob.tol;
        prob = FluenceMapOpt(structs,'x0',prob.x,'tol',prob.tol);   
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