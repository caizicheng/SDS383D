clear classes;close all;clc
numClusters = 100;
numFeatures = 50;
mnistRaw = table2array(readtable('mnist.csv'))';
%mnistRawPlot = mnistRaw;
mnistRaw = mnistRaw - repmat(mean(mnistRaw,2),1,size(mnistRaw,2));
%mnistRaw = mnistRaw ./ repmat(std(mnistRaw')',1,size(mnistRaw,2));
mnistRaw(isnan(mnistRaw)) = 0;
[coeff,score,latent,tsquared,explained,mu] = pca(mnistRaw');
featureSel = coeff(:,1:numFeatures)';
mnistPCA = (featureSel * mnistRaw);
stdMat = eye(numFeatures);%diag(std(mnistPCA'));

mnistPCAReconstrct = inv(coeff') * [mnistPCA;zeros(784-numFeatures,size(mnistRaw,2))];
%mnistPCAReconstrct = mnistPCAReconstrct ./ repmat(std(mnistPCAReconstrct')',1,size(mnistPCAReconstrct,2));
%mnistPCAReconstrct(isnan(mnistPCAReconstrct)) = 0;
figure()

plotSel = 200;
subplot(1,2,1)
imshow(reshape(mnistRaw(:,plotSel),28,28));
title('Original')
subplot(1,2,2)
imshow(reshape(mnistPCAReconstrct(:,plotSel),28,28));
title('PCA Reconstruction Results')

mu = {};
N = size(mnistRaw,2);
mu0 = mean(mnistPCA,2)';
mu0 = repmat(mu0,numClusters,1);
mu{1} = mvnrnd(mu0,stdMat);
Ncount = zeros(1,numClusters);
Ncount(1) = N;
Z = randi(numClusters,N,1);
alpha = 1;
X = mnistPCA';
X = X./repmat(std(X),size(X,1),1);

maxIter = 10;
Mcount = zeros(N,numClusters);
for iter = 1:maxIter

		Mcount = Ncount - (Z(:,end) == 1:numClusters) + alpha;
		for index = 1:N
			for k = 1:numClusters
				
				prob(index,k) = Mcount(index,k)*mvnpdf(X(index,:),mu{end}(k,:),stdMat);
				
			end
			
		end
		
		
		prob = prob./sum(prob,2);
        
		Z(:,end + 1) = zeros(N,1);
		for index = 1:N
			Z(index,end) = sample2(prob(index,:),1);
		end

		XSum = zeros(numClusters,numFeatures);
		NCount = zeros(numClusters,1);
		mu_temp = zeros(numClusters,numFeatures);

		
		for k = 1:numClusters
			XSum(k,:) = sum(X(Z(:,end) == k,:));
			Ncount(k) = sum(Z(:,end) == k);
			s2 = stdMat ./ (Ncount(k) + 1);
			m = XSum(k,:) ./ (Ncount(k) + 1);
			mu_temp(k,:) = mvnrnd(m,s2);			
		end

		mu{end+1} = mu_temp;
		

end
Z_result = Z(:,end);

cluster = {};
clusterSize = zeros(1,numClusters);
for index = 1:numClusters

	cluster{index} = find(Z_result == index);
	clusterSize(index) = length(cluster{index});

end

clusterSel = 9;
figure()
hold on
maxPlot = 4;
for index = 1:maxPlot^2
	subplot(maxPlot,maxPlot,index)
	imshow(reshape(mnistRaw(:,cluster{clusterSel}(index)),28,28)')
end

[clusterSize,clusterIdx] = sort(clusterSize,'descend');
clusterSize = clusterSize(clusterSize > 0);