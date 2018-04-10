clear;
addpath(genpath('kakearney-boundedline-pkg-8179f9a'))
rawData = readtable('faithful.csv','ReadVariableNames',true);
y = rawData.eruptions;
x = rawData.waiting;
numTest = 50;
xstar = linspace(min(x),max(x),numTest)';
N = length(x);
X = [ones(N,1),x,x.^2,x.^3];
Xstar =[ones(numTest,1),xstar,xstar.^2,xstar.^3];
numFeatures = size(X,2);
mu_beta0 = zeros(numFeatures,1);
sigma_beta0 = diag(ones(numFeatures,1));
sigma = 1;
A = sigma^(-2)*X'*X+inv(sigma_beta0);
fstar_mean = 1/sigma^2*Xstar*inv(A)*X'*y;
fstar_covar = Xstar * inv(A) * Xstar';
CI_left = fstar_mean - 1.96*sqrt(diag(fstar_covar));
CI_right = fstar_mean + 1.96*sqrt(diag(fstar_covar));
figure(1)
plot(x,y,'LineStyle','none','Marker','.','MarkerSize',20)
hold on
xlabel('Waiting')
ylabel('Eruptions')
boundedline(xstar, fstar_mean, 1.96*sqrt(diag(fstar_covar)), 'alpha');
