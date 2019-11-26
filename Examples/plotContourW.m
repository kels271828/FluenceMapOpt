function plotContourW(f,wLim,wStep)

figure()

myRed = [0.6350 0.0780 0.1840];
myBlue = [0 0.4470 0.7410];
myOrange = [0.8500 0.3250 0.0980];

% Compute objective function values
wVals = wLim(1):wLim(2);
fVals = zeros(length(wVals));
for i = 1:length(wVals)
    for j = 1:length(wVals)
        f.structs{2}.terms{1}.w = [wVals(i); wVals(j)];
        f.projX();
        fVals(j,i) = 0.5*norm(f.A*f.x - f.d)^2;
    end
end

% Local and global minima
voxels = [1675228; 1674687; 1675607];
load('Gantry16_Couch0_D.mat');
A = D(voxels,86);
load('Gantry352_Couch0_D.mat');
A = [A D(voxels,85)];
lambda = f.lambda;
x1 = (20*lambda*A(2,1) + (81*A(2,2) - 20*A(1,2))*(A(1,1)*A(2,2) -...
    A(1,2)*A(2,1)))/((A(1,1)*A(2,2) - A(1,2)*A(2,1))^2 +...
    lambda*(A(2,1)^2 + A(2,2)^2));
x2 = (20 - A(2,1)*x1)/A(2,2);
wL = f.structs{2}.A*[x1;x2] - 20;
x1 = (20*lambda*A(3,1) + (81*A(3,2) - 20*A(1,2))*(A(1,1)*A(3,2) -...
    A(1,2)*A(3,1)))/((A(1,1)*A(3,2) - A(1,2)*A(3,1))^2 +...
    lambda*(A(3,1)^2 + A(3,2)^2));
x2 = (20 - A(3,1)*x1)/A(3,2);
wG = f.structs{2}.A*[x1;x2] - 20;

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
f.x = f.xInit;
s = strcmp(f.structs{2}.terms{1}.type,'ldvc');
temp = f.structs{2}.A*f.x - f.structs{2}.terms{1}.d;
f.structs{2}.terms{1}.w = (-1)^s*f.projW((-1)^s*temp,f.structs{2}.terms{1}.k);
plot(f.structs{2}.terms{1}.w(1),f.structs{2}.terms{1}.w(2),'o',...
        'MarkerFaceColor',myBlue,'MarkerEdgeColor',myBlue,'MarkerSize',10)
plot(temp(1),temp(2),'s','MarkerFaceColor',myOrange,'MarkerEdgeColor',myOrange)
for i = 1:f.maxIter
    f.projX();
    Axmb = f.structs{2}.A*f.x - f.structs{2}.terms{1}.d;
    coeff = f.structs{2}.terms{1}.step*f.structs{2}.terms{1}.weight/f.structs{2}.nVoxels;
    temp = f.structs{2}.terms{1}.w + coeff*(Axmb - f.structs{2}.terms{1}.w);
    f.structs{2}.terms{1}.w = (-1)^s*f.projW((-1)^s*temp,f.structs{2}.terms{1}.k);
    plot(f.structs{2}.terms{1}.w(1),f.structs{2}.terms{1}.w(2),'o',...
        'MarkerFaceColor',myBlue,'MarkerEdgeColor',myBlue,'MarkerSize',10)
    plot(temp(1),temp(2),'s','MarkerFaceColor',myOrange,'MarkerEdgeColor',myOrange)
end
