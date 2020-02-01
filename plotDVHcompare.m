
% Calculate and plot dose-volume histogram of solution, comparing results from two methods.
function plotDVHcompare(f,xa,xb)

    myLines = lines;

    % Calculate dose-volume histograms
    doses = linspace(0,100,1000);
    dvhInit = zeros(f.nStructs,length(doses));
    dvhFinal_a = zeros(f.nStructs,length(doses));
    dvhFinal_b = zeros(f.nStructs,length(doses));
    for i = 1:f.nStructs
        doseInit = f.structs{i}.A*f.x0; % change back to xInit?
        doseFinal_a = f.structs{i}.A*xa;
        doseFinal_b = f.structs{i}.A*xb;
        for j = 1:length(doses)
            dvhInit(i,j) = 100*sum(doseInit > doses(j))/f.structs{i}.nVoxels;
            dvhFinal_a(i,j) = 100*sum(doseFinal_a > doses(j))/f.structs{i}.nVoxels;
            dvhFinal_b(i,j) = 100*sum(doseFinal_b > doses(j))/f.structs{i}.nVoxels;
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
                if j == 1
                    plot(doses,dvhInit(i,:),'--','Color',[0.4 0.4 0.4])
                    lineHandle_a = plot(doses,dvhFinal_a(i,:),'Color',myLines(2*i-1,:));
                    lineHandle_b = plot(doses,dvhFinal_b(i,:),'Color',myLines(2*i,:));
                    lineName_a = strcat(f.structs{i}.name,' a');
                    lineName_b = strcat(f.structs{i}.name,' b');
                    legendHandles = [legendHandles lineHandle_a lineHandle_b];
                    legendNames = [legendNames, lineName_a lineName_b];
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