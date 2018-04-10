clear;
addpath(genpath('kakearney-boundedline-pkg-8179f9a'))
rawData = readtable('faithful.csv','ReadVariableNames',true);
rawData = rawData(1:10,:);
y = rawData.eruptions;
x = rawData.waiting;
PhiX = [ones(length(x),1),x,x.^2,x.^3];
J = @(ax) -logPFun(PhiX,y,ax(1),ax(2),ax(3),x);
figure(1)
counter = 1;
for row = 1:2
	for col = 1:2

		results = fmincon(J,[abs(rand() * 10) abs(rand() * 10)*1e9 abs(rand() * 10)],[],[],[],[],[0 0 0],[]);

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
		titleStr = sprintf('\\alpha^2 = %0.1e, l^2 = %0.1e, \\sigma^2 = %0.1e',alpha2,l2,sigma2);
		subplot(2,2,counter)
		hold on
		plot(x,y,'LineStyle','none','Marker','.','MarkerSize',20)
		plot(Xstar(1:10:end), fstarbar(1:10:end),'LineStyle','none','Marker','.','MarkerSize',20);
		xlabel('Waiting')
		ylabel('Eruptions')
		boundedline(Xstar, fstarbar, errbound, 'alpha');
		grid on
		title(titleStr)
		counter = counter + 1;
	end
end