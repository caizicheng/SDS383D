clear;
rawData = readtable('restaurants.csv');
X = rawData.Profit;
X = X - mean(X);
X = X / std(X);
[X,idx] = sort(X);
Y = rawData.DinnerService;
Y = Y(idx);
N = length(Y);
Z = randi(2,N,1)-1;
N1 = sum(Z(:,end) == 1);
N0 = length(Y) - N1;

mu0_0 = 0.1;
mu0_1 = 0.9;
sigma0_0 = 1;
sigma0_1 = 1;
sigma2 = 1;
alpha = 1;
beta = 1;
mu0 = [0.5];
mu1 = [0.5];

Prec = X'*X + eye(size(X,2));
maxIter = 100;
for index = 1:maxIter	
	

	%Pz1 = normpdf(X,mu1(end),sigma2(end));
    %Pz0 = normpdf(X,mu0(end),sigma2(end));
    Pz1 = (N1+alpha)/(N-1+alpha+beta)*normpdf(X,mu1(end),sigma2(end));
    Pz0 = (N0+alpha)/(N-1+alpha+beta)*normpdf(X,mu0(end),sigma2(end));

    Pz1 = Pz1./(Pz1 + Pz0);

	r=rand(length(X),1);

	Z(:,end+1)=(r<Pz1);

	N1 = sum(Z(:,end) == 1);
	N0 = length(Y) - N1;

	s2_0 = inv(1/sigma0_0^2+N0/sigma2(end));
	m0 = (mu0(end)/sigma0_0^2+sum(X(Z(:,end) == 0))/sigma2(end))*s2_0;
	mu0(end+1) = normrnd(m0,s2_0);

	s2_1 = inv(1/sigma0_0^2+N1/sigma2(end));
	m1 = (mu1(end)/sigma0_1^2+sum(X(Z(:,end) == 1))/sigma2(end))*s2_1;
	mu1(end+1) = normrnd(m1,s2_1);

	sigma2(end+1)=0.4;

	%p_w = gamrnd(alpha + (size(X,1) + size(X,2)) / 2,(beta+(1/2)*(Y'*Y-mu0(index)*Prec*mu0(index))),1);

end

sum(abs(Z(:,end) - Y))