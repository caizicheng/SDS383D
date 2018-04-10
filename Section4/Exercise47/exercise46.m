clear;
addpath(genpath('kakearney-boundedline-pkg-8179f9a'))
rawData = readtable('faithful.csv','ReadVariableNames',true);
y = rawData.eruptions;
x = rawData.waiting;
PhiX = [ones(length(x),1),x,x.^2,x.^3];
alpha2 = 1;
l2 = 2.5e9;
sigma2 = 1;
alpha2log = alpha2;
l2log = l2;
sigma2log = sigma2;
LRalpha2 = 0.01;
LRl2 = 1e9;
LRsigma2 = 0.01;
K_gen = @(x1,x2) alpha2*exp(-1/2/l2 * norm(x1-x2,2)^2);
for iter = 1:100
	iter	
	K__ = zeros(length(x),length(x));
	partialkpartialalpha2 = K__;
	partialkpartiall2 = K__;
	partialkpartialsigma2 = eye(length(K__));
	for row = 1:length(x)
		for col = 1:length(x)
			K__(row,col) = K_gen(PhiX(row,:),PhiX(col,:)) + sigma2*(row == col);
			partialkpartialalpha2(row,col) = exp(-1/2/l2*norm(PhiX(row,:) - PhiX(col,:),2)^2);
			partialkpartiall2(row,col) = K_gen(PhiX(row,:),PhiX(col,:)) /2/l2^2 * norm(PhiX(row,:) - PhiX(col,:),2)^2;
		end
	end
	partialPpartialalpha2 = 0.5*trace((inv(K__) * y) * (inv(K__) * y)' * partialkpartialalpha2);
	partialPpartiall2 = 0.5*trace((inv(K__) * y) * (inv(K__) * y)' * partialkpartiall2);
	partialPpartialsigma2 = 0.5*trace((inv(K__) * y) * (inv(K__) * y)' * partialkpartialsigma2);
	alpha2 = alpha2 - LRalpha2 * partialPpartialalpha2;
	alpha2log = [alpha2log,alpha2];
	l2 = l2 - LRl2 * partialPpartiall2;
	l2log = [l2log,l2];
	sigma2 = sigma2 - LRsigma2 * partialPpartialsigma2;
	sigma2log = [sigma2log,sigma2];
end

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

fstarbar = Kstar_*inv(K__ + eye(length(K__)) * sigma2) * y;
covfstar = Kstarstar - Kstar_*inv(K__ + eye(length(K__)) * sigma2) * K_star;
errbound = 1.96 * sqrt(diag(covfstar));
titleStr = sprintf('l2 = %0.1f',l2);
figure(1)
plot(x,y,'LineStyle','none','Marker','.','MarkerSize',20)
hold on
xlabel('Waiting')
ylabel('Eruptions')
boundedline(Xstar, fstarbar, errbound, 'alpha');
grid on
title(titleStr)