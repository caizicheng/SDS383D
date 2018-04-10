clear;
addpath(genpath('kakearney-boundedline-pkg-8179f9a'))
rawData = readtable('iris.csv','ReadVariableNames',true);
rawData = rawData(~strcmp(rawData.Species,'setosa'),:);

y = zeros(size(rawData,1),1);
y(strcmp(rawData.Species,'versicolor')) = 1;

x1 = rawData.Sepal_Length;
x2 = rawData.Sepal_Width;
x3 = rawData.Petal_Length;
x4 = rawData.Petal_Width;

PhiX = [x1,x2,x3,x4];

alpha = 1;
l = 5;
sigma = 1;

K_gen = @(xaaa,xbbb) alpha^2*exp(-1/2/l^2 * norm(xaaa-xbbb,2)^2);

PhiXStar = PhiX;

Kstar_ = zeros(length(x1),length(x1));
for row = 1:length(x1)
	for col = 1:length(x1)
		Kstar_(row,col) = K_gen(PhiXStar(row,:),PhiX(col,:));
	end
end
K_star = Kstar_';
Kstarstar = zeros(length(x1),length(x1));
for row = 1:length(x1)
	for col = 1:length(x1)
		Kstarstar(row,col) = K_gen(PhiXStar(row,:),PhiXStar(col,:));
	end
end
K__ = zeros(length(x1),length(x1));
for row = 1:length(x1)
	for col = 1:length(x1)
		K__(row,col) = K_gen(PhiX(row,:),PhiX(col,:));
	end
end

fstarbar = Kstar_*inv(K__ + eye(length(K__)) * sigma^2) * y;

sigma_f = @(f) 1/(1+exp(-f));

logP = @(f) y'*log(sigma_f(f)) - (1-y)'*log(1-sigma_f(f)) - 0.5 *f'*inv(K__)*f - 0.5*log(det(K__));