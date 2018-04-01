%Morgan Kelley
%SDS 383D
%Problem 3.2
clear all
load titanic.mat
n=1312;
%% Data processing
Age=table2array(titanic(2:n+1,3));
Survival=string(table2cell(titanic(2:n+1,5)));
Survival_num=2*ones(n,1);
for i=1:n
    if Survival(i)=='No'
        Survival_num(i)=0;
    else
        Survival_num(i)=1;
    end
end
X=Age;
%Remove values with NA for age
X=Age(isnan(Age)==0);

NA2=0;
    for i=1:n
        if isnan(Age(i))
            NA=i;
            NA2=[NA2; NA];
        end 
    end
  NA2(1)=[];  
    Y=Survival_num;
    Y(NA2)=[]; %Delete corresponding values of Y
 %% Minimize the objective
 lambda=1; t=1; sigma=2;
 fun=@(beta)-((-beta^2/(2*sigma))+sum(-Y.*log(1+exp(-beta.*X))-(1-Y).*log(1+exp(beta.*X)))+lambda*(beta^2-t)); 
 betaR=fminsearch(fun,0); %minimize the negative (maximize)
 %% Plot Problem 3.3
betav=-.1:0.0001:.05;
FF=zeros(size(betav));
for i=1:length(betav)
    FF(i)=fun(betav(i));
end
plot(betav,-FF/max(FF))
%% Mean and precision Problem 3.4
Mean=betaR;
C=-(-(1/sigma^2)+sum(X.^2*exp(sum(X.*betaR))/(1+exp(sum(X.*betaR))).^2.*(2*Y-1)));
%% Problem 3.5 Hessian
%Find Beta_o
X2=[.5*ones(length(X),1) X];
fun2=@(beta)((-(beta*beta')/(2*sigma)) - sum(Y'*log(1+exp(sum(X2*beta'))))-sum((1-Y)'*log(1+exp(sum(X2*beta'))))+lambda*(beta*beta'-t)); 
[betaR2]=fminsearch(fun2,[0,0]);
betav2=[[-.005:.0001:.005]' [-.05:.001:.05]'];
for i=1:length(betav2)
    FF2(i)=fun2(betav2(i,:));
end
%Calculate the Hessian
H11=(1/sigma^2)+sum((X2(:,1).^2.*exp(sum(X2.*betaR2)))/(1+exp(sum(X2.*betaR2))).^2);
H12=sum(Y.*X2.^2.*exp(sum(X2.*betaR2))/(1+exp(sum(X2.*betaR2))).^2);
H21=H12;
H22=(1/sigma^2)+sum((X2(:,2).^2.*exp(sum(X2.*betaR2)))/(1+exp(sum(X2.*betaR2))).^2);
