clear all;
T = readtable('titanic.csv','ReadVariableNames',true);
dataSelect = 1 : size(T,1);
dataSelect = dataSelect(~strcmp(T.Age,'NA'));
data = T(dataSelect,:);
X = zeros(size(data,1),1);
for index = 1:size(data,1)
	X(index) = str2num(data.Age{index});
end

Y = zeros(size(X,1),1);
Y(strcmp(data.Survived,'Yes')) = 1;
%% 3.2
lambda=1; t=1; sigma=2;
p_beta = @(beta) ((-beta^2/(2*sigma))+sum(-Y.*log(1+exp(-beta.*X))-(1-Y).*log(1+exp(beta.*X)))+lambda*(beta^2-t));
obj = @(beta) -p_beta(beta);
beta_MAP = fmincon(obj,0.5,[],[]);
%% 3.3
betaRange = -0.1:0.00001:0.07;
pval = zeros(size(betaRange));
for index = 1:length(betaRange)
	pval(index) = p_beta(betaRange(index));
end
plot(betaRange,pval / abs(max(pval)))
%% 3.5
X2 = [0.5 * ones(size(X,1),1),X];
%p_beta2 = @(beta) -(-1/2/sigma^2 * beta*beta' + sum(X2*beta'.*Y - log(exp(X2*beta') + 1)))% + lambda*(beta * beta' - t));
p_beta2 = @(beta) -((-beta*beta'/(2*sigma))+sum(-Y.*log(1+exp(-X2*beta'))-(1-Y).*log(1+exp(X2*beta')))+lambda*(beta * beta'-t));

beta2_MAP = fmincon(p_beta2,[0 0],[],[])

H = zeros(2,2);
H(1,1) = -1\sigma^2 - sum((X2(:,1).^2 .* exp(X2*beta2_MAP')) ./ (1 + exp(X2 * beta2_MAP')).^2);
H(2,2) = -1\sigma^2 - sum((X2(:,2).^2 .* exp(X2*beta2_MAP')) ./ (1 + exp(X2 * beta2_MAP')).^2);
H(1,2) = -sum((X2(:,1) .* X2(:,2) .* exp(X2*beta2_MAP') ) ./ (1 + exp(X2 * beta2_MAP')).^2);
H(2,1) = H(1,2);