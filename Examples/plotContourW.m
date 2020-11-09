function plotContourW(prob,structs,wLim,wStep)

figure()

myRed = [0.6350 0.0780 0.1840];
myBlue = [0 0.4470 0.7410];
myOrange = [0.8500 0.3250 0.0980];

% Compute objective function values
wVals = wLim(1):wLim(2);
fVals = zeros(length(wVals));
[A,~,~,~] = prob.getA('full');
for ii = 1:length(wVals)
    for jj = 1:length(wVals)
        structs{2}.terms{1}.w = [wVals(ii); wVals(jj)];
        prob.updateStructs(structs,prob.x0);
        d = prob.getd('full');
        prob.x = prob.projX('full');
        fVals(jj,ii) = 0.5*norm(A*prob.x - d)^2;
    end
end

% Local and global minima
voxels = [1675228; 1674687; 1675607];
load('Gantry16_Couch0_D.mat');
A = D(voxels,86);
load('Gantry352_Couch0_D.mat');
A = [A D(voxels,85)];
lambda = prob.lambda;
x1 = (20*lambda*A(2,1) + (81*A(2,2) - 20*A(1,2))*(A(1,1)*A(2,2) -...
    A(1,2)*A(2,1)))/((A(1,1)*A(2,2) - A(1,2)*A(2,1))^2 +...
    lambda*(A(2,1)^2 + A(2,2)^2));
x2 = (20 - A(2,1)*x1)/A(2,2);
wL = prob.structs{2}.A*[x1;x2] - 20;
x1 = (20*lambda*A(3,1) + (81*A(3,2) - 20*A(1,2))*(A(1,1)*A(3,2) -...
    A(1,2)*A(3,1)))/((A(1,1)*A(3,2) - A(1,2)*A(3,1))^2 +...
    lambda*(A(3,1)^2 + A(3,2)^2));
x2 = (20 - A(3,1)*x1)/A(3,2);
wG = prob.structs{2}.A*[x1;x2] - 20;

% Plot objective function contours
contour(wVals,wVals,fVals,15,'LineWidth',2,'LineColor',[0.4 0.4 0.4]), hold on
ax = gca;
ax.XTick = wLim(1):wStep:wLim(2);
ax.YTick = wLim(1):wStep:wLim(2);
ax.XTickLabel = [];
ax.YTickLabel = [];
ax.LineWidth = 2;
axis equal

% Boundary of nonfeasible region
plot([0 wLim(2)],[0 0],'k--','LineWidth',2)
plot([0 0],[0 wLim(2)],'k--','LineWidth',2)

% Plot local and global minima
plot(wL(1),wL(2),'o','MarkerFaceColor',myRed,'MarkerEdgeColor',myRed,'MarkerSize',10)
plot(wG(1),wG(2),'o','MarkerFaceColor',myRed,'MarkerEdgeColor',myRed,'MarkerSize',10)

% Plot iterates of w
prob.x = prob.x0;
k = prob.structs{2}.terms{1}.k;
step = prob.structs{2}.terms{1}.step;
coeff = step*prob.structs{2}.terms{1}.weight/prob.structs{2}.nVoxels;
dose = prob.structs{2}.A*prob.x;
res = dose - prob.structs{2}.terms{1}.d;
structs{2}.terms{1}.w = prob.projW(res,k);
prob.updateStructs(structs,prob.x0);
plot(prob.structs{2}.terms{1}.w(1),prob.structs{2}.terms{1}.w(2),'o',...
        'MarkerFaceColor',myBlue,'MarkerEdgeColor',myBlue,'MarkerSize',10)
plot(res(1),res(2),'s','MarkerFaceColor',myOrange,'MarkerEdgeColor',myOrange)
for ii = 1:prob.maxIter
    prob.x = prob.projX('full');
    dose = prob.structs{2}.A*prob.x;
    res = dose - prob.structs{2}.terms{1}.d;
    wPrev = prob.structs{2}.terms{1}.w;
    wStep = wPrev + coeff*(res - wPrev);
    structs{2}.terms{1}.w = prob.projW(wStep,k);
    prob.updateStructs(structs,prob.x0);
    plot(prob.structs{2}.terms{1}.w(1),prob.structs{2}.terms{1}.w(2),'o',...
        'MarkerFaceColor',myBlue,'MarkerEdgeColor',myBlue,'MarkerSize',10)
    plot(wStep(1),wStep(2),'s','MarkerFaceColor',myOrange,'MarkerEdgeColor',myOrange)
end
