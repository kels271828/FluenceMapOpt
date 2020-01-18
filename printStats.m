
function printStats(f,x)
    
    for i = 1:f.nStructs
        fprintf('Structure: %s\n',f.structs{i}.name)
        for j = 1:length(f.structs{i}.terms)
            if strcmp(f.structs{i}.terms{j}.type,'unif') % Uniform PTV objective
                dose = f.structs{i}.A*x;
                obj = norm(dose - f.structs{i}.terms{j}.d);
                D95 = percentile(dose,0.95);
                fprintf('* unif | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f\n',...
                    min(dose), mean(dose), median(dose), D95, max(dose), std(dose), obj);  
            elseif strcmp(f.structs{i}.terms{j}.type,'ldvc') % LDVC PTV objective
                fprintf('* ldvc | %.2f\n', 100-getPercent(f,i,j,x));
            else % UDVC OAR objective
                fprintf('* udvc | %.2f | %.2f\n', getPercent(f,i,j,x), dvhArea(f,i,x))
            end
        
        end
    end

end