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


K = zeros(length(x1),length(x1));
for row = 1:length(x1)
	for col = 1:length(x1)
		K(row,col) = K_gen(PhiX(row,:),PhiX(col,:));
	end
end

K = K + sigma^2*eye(length(y));
invK = inv(K);
obj = @(f) objFun(y,f,K);
options = optimoptions('fmincon','Algorithm','interior-point','UseParallel',true,'Display','iter-detailed');
fMAP = fmincon(obj,randn(length(y),1),[],[],[],[],[],[],[],options);

H = zeros(length(y));
for index = 1:length(y)
	p_y_f = ((1/(1+exp(-fMAP(index))))^y(index))*((1-1/(1+exp(-fMAP(index))))^(1-y(index)));
	partialPpartialf = y(index)*exp(fMAP(index))/(1+exp(fMAP(index))) - exp(fMAP(index)*(y(index)+1))/(1+exp(fMAP(index)))^2 - fMAP(index)*sum(invK(index,:));
	partialP2partialf2 = y(index)^2*exp(fMAP(index))/(exp(fMAP(index))+1)-exp(fMAP(index))^(y(index)+1)/(exp(fMAP(index))+1)^2+2*exp(fMAP(index))^(y(index)+2)/(exp(fMAP(index))+1)^3;
	H(index,index) = -1/p_y_f^2*partialPpartialf^2+1/p_y_f*partialP2partialf2-sum(invK(index,:));
end
boundErr = sqrt(abs(diag(-inv(H))));

p = 1./(1+exp(-fMAP));
pl = 1./(1+exp(-fMAP + 1.96 * boundErr));
ph = 1./(1+exp(-fMAP - 1.96 * boundErr));
figure(1)
hold on
plot(p)
plot(y)
plot(pl)
plot(ph)
legend('Probability','Truth','High Bound','Low Bound')