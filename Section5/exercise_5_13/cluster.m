function [ID] = cluster(alpha,n)
id=zeros(n,1);
K=10; %clusters
for k=1:5
    P=drchrnd(alpha*ones(n,1)');
    keyboard
    id=id+sum(mnrnd(1,P,K),1)';
end
ID=find(id>0);
end

