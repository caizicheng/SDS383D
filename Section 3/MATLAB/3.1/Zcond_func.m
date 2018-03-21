function [Z,RHO,NAR,NGIBBS] = Zcond_func(beta,tau,X,Y)
%Conditional z generator based on truncated multivariate normal
%
mean=(X*beta');
stddev=sqrt(1/tau)*eye(size(X,1),size(X,1));
A=zeros(1,768);
for i=1:length(Y)
if Y(i)==1
   A(1,i)=1;
elseif Y(i)==0
   A(1,i)=-1;
end
b=0;
end
[Z,RHO,NAR,NGIBBS]=rmvnrnd_func(mean,stddev,1,A,b);
end

