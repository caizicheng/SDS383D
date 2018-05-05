%Morgan Kelley
%Ex. 5.7
clear all; clc;
%%
Ns=10;
N=100;
K=5;

z_true = [ones(1,N) 2*ones(1,N) 3*ones(1,N) 4*ones(1,N) 5*ones(1,N)];
N2 = length(z_true);
X=zeros(N2,2);
mu_true = ones(5,2);
mu_true(2,:) = 2;
mu_true(3,:) = 3;
mu_true(4,:) = 4;
mu_true(5,:) = 5;

mu_true = 4*mu_true;
for i=1:N2
    X(i,:)=mvnrnd(mu_true(z_true(i),:),eye(2));
end

mu0=mu_true; sigma0=std(X);

[mut,zt,Ncountt]=Ex57Func_7(Ns,X,K,mu0,sigma0);
%%
%%
clf
figure(1); hold on
for i=1:N
plot(X(i,1),X(i,2),'g.','MarkerSize',8)
end

for i=i+1:N*2
plot(X(i,1),X(i,2),'k.','MarkerSize',8)
end

for i=i+1:N*3
plot(X(i,1),X(i,2),'r.','MarkerSize',8)
end

for i=i+1:N*4
plot(X(i,1),X(i,2),'b.','MarkerSize',8)
end

for i=i+1:N2
plot(X(i,1),X(i,2),'m.','MarkerSize',8)
end

figure(2)
plot(zt(end,:),'.');

figure(3)
histogram(zt(end,:));

Error=sum(abs(z_true-zt(Ns+1,:)));
PE=Error/N2;
