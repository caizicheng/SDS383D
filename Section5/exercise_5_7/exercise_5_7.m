clear classes;close all;clc
numClusters = 5;
mu_true = linspace(1,50,numClusters)';
samples = {};
nSamples = 100;
N = nSamples * numClusters;
mu0 = [mu_true mu_true];
mu{1} = mvnrnd(mu0,eye(size(mu0,2)));
Ncount = zeros(1,numClusters);
Ncount(1) = N;
Z = randi(5,N,1);
alpha = 1;
X = [];
%figure()
%hold on
for index = 1:numClusters
	samples{index} = randn(nSamples,2) + mu_true(index);
	X = [X;samples{index}];
	%plot(samples{index}(:,1),samples{index}(:,2),'LineStyle','None','Marker','.')
end
%hold off

maxIter = 10;
Mcount = zeros(N,numClusters);
for iter = 1:maxIter

		Mcount = Ncount - (Z(:,end) == 1:numClusters) + alpha;
		for index = 1:N
			for k = 1:numClusters
				prob(index,k) = Mcount(index,k)*mvnpdf(X(index,:),mu{end}(k,:),eye(size(mu0,2)));
			end
		end
		

		prob = prob./sum(prob,2);
		Z(:,end + 1) = zeros(N,1);
		for index = 1:N
			Z(index,end) = sample2(prob(index,:),1);
		end
		XSum = zeros(numClusters,2);
		NCount = zeros(numClusters,1);
		mu_temp = zeros(numClusters,2);

		for k = 1:numClusters
			XSum(k,:) = sum(X(Z(:,end) == k));
			Ncount(k) = sum(Z(:,end) == k);
			s2 = eye(2) ./ (Ncount(k) + 1);
			m = XSum(k,:) ./ (Ncount(k) + 1);
			mu_temp(k,:) = mvnrnd( m,s2 );			
		end
		mu{end+1} = mu_temp;
		

end
Z_result = Z(:,end);
results = {};
figure()
hold on
for index = 1:numClusters
	results{index} = X(Z_result == index,:);
	plot(results{index}(:,1),results{index}(:,2),'LineStyle','None','Marker','.');
end
hold off
