clear;
addpath(genpath('kakearney-boundedline-pkg-8179f9a'))
rawData = readtable('faithful.csv','ReadVariableNames',true);
y = rawData.eruptions;
x = rawData.waiting;
PhiX = [ones(length(x),1),x,x.^2,x.^3];

J = @(ax) -logPFun(PhiX,y,ax(1),ax(2),ax(3),x);
results = fminsearch(J,[1 2.5e9 1]);

alpha2 = results(1);
l2 = results(2);
sigma2 = results(3);
K_gen = @(x1,x2) alpha2*exp(-1/2/l2 * norm(x1-x2,2)^2);
Xstar = linspace(0,100,200)';
PhiXStar = [ones(length(Xstar),1),Xstar,Xstar.^2,Xstar.^3];
Kstar_ = zeros(length(Xstar),length(x));
for row = 1:length(Xstar)
	for col = 1:length(x)
		Kstar_(row,col) = K_gen(PhiXStar(row,:),PhiX(col,:));
	end
end
K_star = Kstar_';
Kstarstar = zeros(length(Xstar),length(Xstar));
for row = 1:length(Xstar)
	for col = 1:length(Xstar)
		Kstarstar(row,col) = K_gen(PhiXStar(row,:),PhiXStar(col,:));
	end
end
K__ = zeros(length(x),length(x));
for row = 1:length(x)
	for col = 1:length(x)
		K__(row,col) = K_gen(PhiX(row,:),PhiX(col,:));
	end
end

fstarbar = Kstar_*inv(K__ + eye(length(K__)) * sigma2) * y;
covfstar = Kstarstar - Kstar_*inv(K__ + eye(length(K__)) * sigma2) * K_star;
errbound = 1.96 * sqrt(diag(covfstar));
titleStr = sprintf('l2 = %0.1f',l2);
figure(1)
hold on
plot(x,y,'LineStyle','none','Marker','.','MarkerSize',20)
plot(Xstar(1:10:end), fstarbar(1:10:end),'LineStyle','none','Marker','.','MarkerSize',20);
xlabel('Waiting')
ylabel('Eruptions')
boundedline(Xstar, fstarbar, errbound, 'alpha');
grid on
title(titleStr)