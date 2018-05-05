%Morgan Kelley
%Ex 5.17
close all; clear all; clc;
%% Load data
load mnist.mat;

x1=mnist-mean(mnist,2);
x2=mnist./std(x1);
x3=x2-mean(x2,2);

%% PCA to reduce dimensionality
elem=50;
Xred=pca(x3','NumComponents',elem);

[eig,scores,lat,ts,ex,M1]=pca(x3','NumComponents',elem);
rec=scores*eig'+repmat(M1,784,1);
rec=rec';
%% Get truths
Ns=10;
K=100;
[N,D]=size(Xred);
N2=N/K;

z_true = [zeros(1,N2) ones(1,N2) 2*ones(1,N2) 3*ones(1,N2) 4*ones(1,N2) 5*ones(1,N2) 6*ones(1,N2) 7*ones(1,N2) 8*ones(1,N2) 9*ones(1,N2)];

%% Run Gibbs Sampler
[mut,zt,Ncountt]=Ex517Func(Ns,Xred,K);

%% Plot
clf
for i=2:Ns+1
    figure(1)
    subplot(2,5,i-1)
histogram(zt(i,:)-1,50);
title(['Iteration: ' num2str(i)])
end
