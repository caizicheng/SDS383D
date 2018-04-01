clear classes;
T = readtable('pima.csv','ReadVariableNames',true);
intercept = ones(size(T,1),1);
X = table2array(T);
%X = [X(:,1:end-1),intercept];
X = X(:,1:end-1);
Y = T.class_variable;

beta = gibbsSampler(X,Y,1000);

betaSelect = mean(beta(:,185:end),2);

z = X * betaSelect;

Yhat = zeros(size(z));
Yhat(z > 0) = 1;
correctRate = mean(Yhat == Y);