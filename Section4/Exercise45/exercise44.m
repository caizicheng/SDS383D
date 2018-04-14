clear;
addpath(genpath('kakearney-boundedline-pkg-8179f9a'))
rawData = readtable('faithful.csv','ReadVariableNames',true);
y = rawData.eruptions;
x = rawData.waiting;
PhiX = [ones(length(x),1),x,x.^2,x.^3];
alpha = 1;
l = 50000;
sigma = 1;
K_gen = @(x1,x2) alpha^2*exp(-1/2/l^2 * norm(x1-x2,2)^2);
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

fstarbar = Kstar_*inv(K__ + eye(length(K__)) * sigma^2) * y;
covfstar = Kstarstar - Kstar_*inv(K__ + eye(length(K__)) * sigma^2) * K_star;
errbound = 1.96 * sqrt(diag(covfstar));
titleStr = sprintf('l = %0.1f',l);
figure(1)
plot(x,y,'LineStyle','none','Marker','.','MarkerSize',20)
hold on
xlabel('Waiting')
ylabel('Eruptions')
boundedline(Xstar, fstarbar, errbound, 'alpha');
grid on
title(titleStr)