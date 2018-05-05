function [mu,z,Ncount] = Ex517Func(Ns,X,K)
[N,D] = size(X);
n=N/K;
%Initial values
mu=cell(Ns+1);
z=ones(Ns+1,N);
z=[];
for k=1:K
    z2(k,:)=k*ones(1,n);
    z=[z z2(k,:)];
end

alpha=1/K;
Ncount=zeros(1,K);
Xsum=zeros(D,K);
for i=1:K
    Ncount(i)=sum(z(1,:)==i);
    Xs(:,i)=sum(X(z(1,:)==i,:),1)';
end
for k=1:K
        s20=eye(D)./(Ncount(1,k)+1);
        m0=Xsum(:,k)./(Ncount(1,k)+1);
        mu{1}(k,:)=mvnrnd( m0,s20 );
    end
%Begin Sampling
for iter=2:Ns+1
    prob=zeros(N,K);
    Mcount=zeros(N,K);
    for i=1:N
        for k=1:K
            Mcount(i,k)=Ncount(k) -double(z(iter-1,i)==k)+ alpha;
            prob(i,k) = Mcount(i,k)*mvnpdf(X(i,:),mu{iter-1}(k,:),eye(D));
        end
        prob(i,:) = prob(i,:)/sum(prob(i,:));
        z(iter,i)= sample2(prob(i,:),1); 
    end
    mu{iter} = zeros(K,D);
    for k=1:K
        Ncount(k)=sum(z(1,:)==k);
        Xs(:,k)=sum(X(z(1,:)==k,:),1)';
        s2=eye(D)./(Ncount(k)+1);%( (1./sigma0.^2) + (D*(Ncount(k))./sigma2) ).^(-1);
        m=Xsum(:,k)./(Ncount(k)+1);%( (mu0(k,:)./sigma0.^2) + ((Xsum(k,:))./sigma2) )* 1./( (1./sigma0.^2) + (D*(Ncount(k))./sigma2) );
        mu{iter}(k,:)=mvnrnd( m,s2 );
    end
end
end

