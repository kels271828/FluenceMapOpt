function plotContourX(lambda,alpha,w,x)

figure()

% Get voxels and dose matrix
voxels = [1675228; 1674687; 1675607];
load('Gantry16_Couch0_D.mat');
A = D(voxels,86);
load('Gantry352_Couch0_D.mat');
A = [A D(voxels,85)];

% Plot contours of objective function
xVals = 0:1e2:2e4;
[X1,X2] = meshgrid(xVals);
Z = 1/2*(A(1,1)*X1 + A(1,2)*X2 - 81).^2 + lambda/2*(X1.^2 + X2.^2);
if nargin > 1
    W1 = w(1) - (A(2,1)*X1 + A(2,2)*X2 - 20);
    W2 = w(2) - (A(3,1)*X1 + A(3,2)*X2 - 20); 
    Z = Z + alpha/4*(W1.^2 + W2.^2);
end
contour(X1,X2,Z,15,'Color',[0.4 0.4 0.4],'LineWidth',2), hold on
axis equal

% Boundary of nonfeasible region
bVals1 = [20*(A(3,2) - A(2,2))/(A(2,1)*A(3,2) - A(3,1)*A(2,2)),...
    (20 - 2.5e4*A(2,2))/A(2,1), 2.5e4, 2.5e4];
bVals2 = [(20 - A(2,1)*(20*(A(3,2) - A(2,2))/(A(2,1)*A(3,2) -...
    A(3,1)*A(2,2))))/A(2,2), 2.5e4, 2.5e4, (20 - 2.5e4*A(3,1))/A(3,2)];
if nargin > 1
    fill(bVals1,bVals2,'w','LineWidth',2,'LineStyle','--','FaceAlpha',0)
else
    fill(bVals1,bVals2,'w','LineWidth',2)
end

% Reset axes limits
ax = gca;
ax.XLim = [min(xVals) max(xVals)];
ax.YLim = [min(xVals) max(xVals)];

% Increase line width
ax.LineWidth = 2;

% Plot local minima
markerColor1 = [0.6350 0.0780 0.1840];
A1 = (20*lambda*A(2,1) + (81*A(2,2) - 20*A(1,2))*(A(1,1)*A(2,2) -...
    A(1,2)*A(2,1)))/((A(1,1)*A(2,2) - A(1,2)*A(2,1))^2 +...
    lambda*(A(2,1)^2 + A(2,2)^2));
A2 = (20 - A(2,1)*A1)/A(2,2);
plot(A1,A2,'o','MarkerFaceColor',markerColor1,...
    'MarkerEdgeColor',markerColor1,'MarkerSize',10)
B1 = (20*lambda*A(3,1) + (81*A(3,2) - 20*A(1,2))*(A(1,1)*A(3,2) -...
    A(1,2)*A(3,1)))/((A(1,1)*A(3,2) - A(1,2)*A(3,1))^2 +...
    lambda*(A(3,1)^2 + A(3,2)^2));
B2 = (20 - A(3,1)*B1)/A(3,2);
plot(B1,B2,'o','MarkerFaceColor',markerColor1,...
    'MarkerEdgeColor',markerColor1,'MarkerSize',10)

% Plot global minimum of current iterate of relaxed problem
if nargin > 1
    markerColor2 = [0 0.4470 0.7410];
    plot(x(1),x(2),'go','MarkerFaceColor',markerColor2,...
        'MarkerEdgeColor',markerColor2,'MarkerSize',10)
end

% Set axes ticks and remove labels
ax.XTick = (0:0.5:2)*1e4;
ax.YTick = (0:0.5:2)*1e4;
ax.XTickLabel = [];
ax.YTickLabel = [];
