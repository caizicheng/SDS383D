clear;
addpath(genpath('kakearney-boundedline-pkg-8179f9a'))
rawData = readtable('iris.csv','ReadVariableNames',true);

y_setosa = zeros(size(rawData,1),1);
y_setosa(strcmp(rawData.Species,'setosa')) = 1;

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

K = K + sigma^2*eye(length(y_setosa));
invK = inv(K);
obj = @(f) objFun(y_setosa,f,K);
options = optimoptions('fmincon','Algorithm','interior-point','UseParallel',true,'Display','iter-detailed');
fMAP_setosa = fmincon(obj,randn(length(y_setosa),1),[],[],[],[],[],[],[],options);

y_versicolor = zeros(size(rawData,1),1);
y_versicolor(strcmp(rawData.Species,'versicolor')) = 1;
obj = @(f) objFun(y_versicolor,f,K);
fMAP_versicolor = fmincon(obj,randn(length(y_versicolor),1),[],[],[],[],[],[],[],options);

y_virginica = zeros(size(rawData,1),1);
y_virginica(strcmp(rawData.Species,'virginica')) = 1;
obj = @(f) objFun(y_virginica,f,K);
fMAP_virginica = fmincon(obj,randn(length(y_virginica),1),[],[],[],[],[],[],[],options);

fMAP = [fMAP_setosa,fMAP_versicolor,fMAP_virginica];
[~,ident] = max(fMAP');
truth = zeros(size(ident));
truth(strcmp(rawData.Species,'setosa')) = 1;
truth(strcmp(rawData.Species,'versicolor')) = 2;
truth(strcmp(rawData.Species,'virginica')) = 3;
plot(ident,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
plot(truth,'LineStyle','none','Marker','o','MarkerSize',10)
figure()
hold on
plot(1./exp(-fMAP_setosa));
plot(1./exp(-fMAP_versicolor));
plot(1./exp(-fMAP_virginica));
legend('setosa','versicolor','virginica')