clear;close all;
alpha = 1;
x = linspace(0,100,200);
l = 10;
K = zeros(length(x),length(x));
for row = 1:length(x)
	for col = 1:length(x)
		K(row,col) = alpha^2 * exp(-1/2/l^2 * norm(x(row) - x(col),2)^2);
	end
end

numSamples = 5;
mu = zeros(length(x),1);
x_samp = zeros(length(x),numSamples);
figure()
hold on
for sampleIndex = 1:numSamples
	x_samp(:,sampleIndex) = mvnrnd(mu,K);
	plot(x_samp(:,sampleIndex))
end
titlestr = sprintf('Gaussian Process l = %0.1f',l);
title(titlestr)