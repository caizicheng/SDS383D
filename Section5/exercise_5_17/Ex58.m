%Morgan Kelley
%Ex 5.8
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
K=10;
[N,D]=size(Xred);
N2=N/K;

z_true = [zeros(1,N2) ones(1,N2) 2*ones(1,N2) 3*ones(1,N2) 4*ones(1,N2) 5*ones(1,N2) 6*ones(1,N2) 7*ones(1,N2) 8*ones(1,N2) 9*ones(1,N2)];

%% Run Gibbs Sampler
[mut,zt,Ncountt]=Ex58Func_2(Ns,Xred,K);

%% Plot
close all;
figure(1)
subplot(1,2,1)
imshow(reshape(mnist(801,:),28,28));
title('Original')
subplot(1,2,2)
imshow(reshape(rec(801,:),28,28));
title('Reconstructed')

figure(2)
plot(z_true,'.'); hold on;
plot(zt(end,:)-1,'.');
ylabel('Digit')
xlabel('Sample')
legend('True','Predicted')

figure(3)
histogram(zt(end,:)-1);
count=0;
for i=1:N
if (z_true(i)==zt(end,i)-1)
    count=1+count;
end
end

Error=1-count/N;
