function [prob,nIter,ttime] = calcBeamsContinue(prob,structs,gamma,...
    sigma,stopType,maxIter,print,consolidate)
    % CALCBEAMSCONTINUE Approach inspired by Lu paper.
    
    % Initialization
    t1 = clock;
    prob.initProb(false);
    obj = prob.getObj();
    d95 = prob.getPercentile(1,0.95);
    structsOld = structs;
    dPos = getPos(prob,structsOld);
    if print
        fprintf('iter: 0, obj: %7.4e, dPos: %4.2f\n',obj,dPos(1,1));
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
            dPos = [dPos; getPos(prob,structsOld)];
        end
            
        % Print status
        if print
            fprintf('iter: %d, obj: %7.4e, dPos: %4.2f\n',kk,obj(end),dPos(end,1));
        end
        
        % Check convergence
        if checkConv(prob,stopType,structsOld,d95)
            break
        end
        
        % Update parameters
        for ii = 1:prob.nStructs
            for jj = 1:prob.structs{ii}.nTerms
                if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                    weight = prob.structs{ii}.terms{jj}.weight;
                    percent = prob.structs{ii}.terms{jj}.percent;
                    dose = prob.structs{ii}.terms{jj}.dose;
                    structs{ii}.terms{jj}.weight = (1 + sigma)*weight;
                    structs{ii}.terms{jj}.percent = (1 - sigma)*percent;
                    if strcmp(prob.structs{ii}.terms{jj}.type,'udvc')
                        structs{ii}.terms{jj}.dose = (1 - sigma)*dose;
                    else
                        structs{ii}.terms{jj}.dose = (1 + sigma)*dose;
                    end
                    prob.updateStructs(structs,prob.x)
                end
            end
        end
        prob.tol = gamma*prob.tol;
    end
    nIter = kk;
    t2 = clock;
    ttime = etime(t2,t1);
end

function dPos = getPos(prob,structs)
    % Get percent voxels above dose value.
    dPos = [];
    for ii = 1:prob.nStructs
        for jj = 1:prob.structs{ii}.nTerms
            if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                temp = prob.getPercent(ii,jj,structs{ii}.terms{jj}.dose);
                dPos = [dPos temp];
            end
        end
    end
end

function stop = checkConv(prob,stopType,structs,d95)
    % Check convergence (meet OAR DVC or violate PTV D95)
    stop = 1;
    if stopType == 0
        for ii = 1:prob.nStructs
            for jj = 1:prob.structs{ii}.nTerms
                if ~strcmp(prob.structs{ii}.terms{jj}.type,'unif')
                    percent = structs{ii}.terms{jj}.percent;
                    temp = prob.getPercent(ii,jj,structs{ii}.terms{jj}.dose);
                    if strcmp(prob.structs{ii}.terms{jj}.type,'udvc')
                        if temp > percent
                            stop = 0;
                        end
                    else
                        if temp < percent
                            stop = 0;
                        end
                    end
                end
            end
        end
    else
        if prob.getPercentile(1,0.95) > 0.98*d95
            stop = 0;
        end
    end
end
