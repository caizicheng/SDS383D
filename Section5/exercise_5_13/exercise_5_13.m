clear all;close all;clc;
alpha = [0.001 0.1 1 10];
numSamples = 100;
for alphaInd = 1:length(alpha)
	id{alphaInd} = zeros(numSamples,1);
	for k = 1:5
		pi = drchrnd(alpha(alphaInd) * ones(1,numSamples));
		id{alphaInd} = id{alphaInd} + sum(mnrnd(1,pi,10))';
	end
end

resultCount = [];
for alphaInd = 1:length(alpha)
	result = find(id{alphaInd}>0);
	resultCount(end+1) = length(result);
end