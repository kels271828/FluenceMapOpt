function prob = calcBeamsContinue(prob,structs,gamma,sigma,maxIter,...
    print,stopType,consolidate)
    % CALCBEAMSTEMP Approach inspired by Lu paper.
    
    % Initialization
    t1 = clock;
    structsOld = structs;
    prob.initProb(false);
    obj = [prob.getObj()];
    d95 = getD95(prob);
    dvc = getDVC(prob,structsOld);
    if print
        printResults(0,obj,d95,dvc);
    end
    
    % Continuation
    for kk = 1:maxIter
        % Solve approximate solution
        if consolidate
            prob = calcBeamsConsolidate(prob,false);
        else
            prob.calcBeams(false);
        end
        obj = [obj, prob.obj(end)];
        d95 = [d95; getD95(prob)];
        dvc = [dvc; getDVC(prob,structsOld)];
        if print
            printResults(kk,obj,d95,dvc);
        end
        
        % Check convergence
        if converged(prob,structsOld,stopType,d95)
            t2 = clock;
            results.nIter = kk;
            results.obj = obj;
            results.time = etime(t2,t1);
            prob.updateResults(results);
            return
        end
        
        % Update parameters
        structs = updatePars(prob,structs,structsOld,sigma,stopType);
        prob.updateStructs(structs,prob.x);
        prob.tol = gamma*prob.tol;
    end
end
    
function d95 = getD95(prob)
    d95 = [];
    for ii = 1:prob.nStructs
        for jj = 1:prob.structs{ii}.nTerms
            if strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                d95 = [d95 prob.getPercentile(ii,0.95)];
            end
        end
    end
end
    
function dvc = getDVC(prob,structs)
    dvc = [];
    for ii = 1:prob.nStructs
        for jj = 1:prob.structs{ii}.nTerms
            if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                dose = structs{ii}.terms{jj}.dose;
                dvc = [dvc prob.getPercent(ii,jj,dose)];
            end
        end
    end
end

function printResults(iter,obj,d95,dvc)
    fprintf('iter: %d, obj: %7.4e\n',iter,obj(end))
    fprintf('d95:')
    fprintf('%6.2f',d95(end,:))
    fprintf('\ndvc:')
    fprintf('%6.2f',dvc(end,:))
    fprintf('\n\n')
end

function stop = converged(prob,structs,stopType,d95)
    if stopType == 0
        stop = true;
        for ii = 1:prob.nStructs
            for jj = 1:prob.structs{ii}.nTerms
                if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                    dose = structs{ii}.terms{jj}.dose;
                    percent = structs{ii}.terms{jj}.percent;
                    dvc = prob.getPercent(ii,jj,dose);
                    if dvc > percent
                        stop = false;
                        return
                    end
                end
            end
        end
    else
        stop = false;
        count = 1;
        for ii = 1:prob.nStructs
            for jj = 1:prob.structs{ii}.nTerms
                if strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                    if prob.getPercentile(ii,0.95) <= 0.98*d95(1,count)
                        stop = true;
                        return
                    end
                end
            end
            count = count + 1;
        end
    end
end

function structs = updatePars(prob,structs,structsOld,sigma,stopType)
    for ii = 1:prob.nStructs
        for jj = 1:prob.structs{ii}.nTerms
            if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                weight = prob.structs{ii}.terms{jj}.weight;
                percent = prob.structs{ii}.terms{jj}.percent;
                dose = prob.structs{ii}.terms{jj}.dose;
                dvc = prob.getPercent(ii,jj,structsOld{ii}.terms{jj}.dose);
                if dvc > structsOld{ii}.terms{jj}.percent || stopType == 1
                    structs{ii}.terms{jj}.weight = (1 + sigma)*weight;
                    structs{ii}.terms{jj}.percent = (1 - sigma)*percent;
                    if strcmp(prob.structs{ii}.terms{jj}.type,'udvc')
                        structs{ii}.terms{jj}.dose = (1 - sigma)*dose;
                    else
                        structs{ii}.terms{jj}.dose = (1 + sigma)*dose;
                    end
                end
            end
        end
    end
end
