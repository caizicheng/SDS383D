clear;
addpath(genpath('kakearney-boundedline-pkg-8179f9a'))
rawData = readtable('weather.csv','ReadVariableNames',true);

y = rawData.temperature;
x1 = rawData.lon;
x2 = rawData.lat;

PhiX = [ones(length(x1),1),x1,x1.^2,x1.^3,x2,x2.^2,x2.^3];
J = @(ax) -logPFun(PhiX,y,ax(1),ax(2),ax(3),x1);

results = fmincon(J,[abs(rand() * 10) abs(rand() * 10)*1e9 abs(rand() * 10)],[],[],[],[],[0 0 0],[]);

alpha2 = results(1);
l2 = results(2);
sigma2 = results(3);
K_gen = @(x1,x2) alpha2*exp(-1/2/l2 * norm(x1-x2,2)^2);
X1star = linspace(-140,-110,40)';
X2star = linspace(35,55,40)';
%PhiXStar = [ones(length(X1star),1),X1star,X1star.^2,X1star.^3,X2star,X2star.^2,X2star.^3];
PhiXStar = [];
for x1index = 1:length(X1star)
	for x2index = 1:length(X2star)
		current = [1,X1star(x1index),X1star(x1index)^2,X1star(x1index)^3,X2star(x2index),X2star(x2index)^2,X2star(x2index)^3];
		PhiXStar = [PhiXStar;current];
	end
end

Kstar_ = zeros(length(PhiXStar),length(x1));
for row = 1:length(PhiXStar)
	for col = 1:length(x1)
		Kstar_(row,col) = K_gen(PhiXStar(row,:),PhiX(col,:));
	end
end
K_star = Kstar_';
Kstarstar = zeros(length(PhiXStar),length(PhiXStar));
for row = 1:length(PhiXStar)
	for col = 1:length(PhiXStar)
		Kstarstar(row,col) = K_gen(PhiXStar(row,:),PhiXStar(col,:));
	end
end
K__ = zeros(length(x1),length(x1));
for row = 1:length(x1)
	for col = 1:length(x1)
		K__(row,col) = K_gen(PhiX(row,:),PhiX(col,:));
	end
end

fstarbar = Kstar_*inv(K__ + eye(length(K__)) * sigma2) * y;
covfstar = Kstarstar - Kstar_*inv(K__ + eye(length(K__)) * sigma2) * K_star;
errbound = 1.96 * sqrt(diag(covfstar));
titleStr = sprintf('\\alpha^2 = %0.1e, l^2 = %0.1e, \\sigma^2 = %0.1e',alpha2,l2,sigma2);

temperaturePost = reshape(fstarbar,length(X1star),length(X2star));
[plotX,plotY] = meshgrid(X1star,X2star);
plotLen = 10;
selectIndex = randi(length(x1),plotLen,1);
figure()
contourf(plotX,plotY,temperaturePost,'ShowText','on')
hold on
scatter(x1(selectIndex),x2(selectIndex))
a = y(selectIndex); b = num2str(a); c = cellstr(b);
dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
text(x1(selectIndex)+dx, x2(selectIndex)+dy, c);
%figure(1)
%hold on
%plot(x,y,'LineStyle','none','Marker','.','MarkerSize',20)
%plot(Xstar(1:10:end), fstarbar(1:10:end),'LineStyle','none','Marker','.','MarkerSize',20);
%xlabel('Waiting')
%ylabel('Eruptions')
%boundedline(Xstar, fstarbar, errbound, 'alpha');
%grid on
%title(titleStr)
