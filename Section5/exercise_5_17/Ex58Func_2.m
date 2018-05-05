function [mu,z,Ncount] = Ex58Func_2(Ns,X,K)
[N,D] = size(X);
n=N/K;
%Initial values
mu=cell(Ns+1);
z=ones(Ns+1,N);
z(1,:)= [2*ones(1,n) 1*ones(1,n) 5*ones(1,n) 4*ones(1,n) 3*ones(1,n) 7*ones(1,n) 9*ones(1,n) 6*ones(1,n) 8*ones(1,n) 10*ones(1,n) ];
alpha=1;
Ncount=zeros(1,K);
Ncount(1,:) = [sum(z(1,:)==1) sum(z(1,:)==2) sum(z(1,:)==3) sum(z(1,:)==4) sum(z(1,:)==5) sum(z(1,:)==6) sum(z(1,:)==7) sum(z(1,:)==8) sum(z(1,:)==9) sum(z(1,:)==10)];
Xsum = [sum(X(z(1,:)==1,:),1)' sum(X(z(1,:)==2,:),1)' sum(X(z(1,:)==3,:),1)' sum(X(z(1,:)==4,:),1)' sum(X(z(1,:)==5,:),1)' sum(X(z(1,:)==6,:),1)' sum(X(z(1,:)==7,:),1)' sum(X(z(1,:)==8,:),1)' sum(X(z(1,:)==9,:),1)' sum(X(z(1,:)==10,:),1)'];
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
    Xsum = [sum(X(z(iter,:)==1,:),1)' sum(X(z(iter,:)==2,:),1)' sum(X(z(iter,:)==3,:),1)' sum(X(z(iter,:)==4,:),1)' sum(X(z(iter,:)==5,:),1)' sum(X(z(iter,:)==6,:),1)' sum(X(z(iter,:)==7,:),1)' sum(X(z(iter,:)==8,:),1)' sum(X(z(iter,:)==9,:),1)' sum(X(z(iter,:)==10,:),1)'];
    Ncount = [sum(z(iter,:)==1) sum(z(iter,:)==2) sum(z(iter,:)==3) sum(z(iter,:)==4) sum(z(iter,:)==5) sum(z(iter,:)==6) sum(z(iter,:)==7) sum(z(iter,:)==8) sum(z(iter,:)==9) sum(z(iter,:)==10)];    
    
    mu{iter} = zeros(K,D);
    for k=1:K
        s2=eye(D)./(Ncount(k)+1);%( (1./sigma0.^2) + (D*(Ncount(k))./sigma2) ).^(-1);
        m=Xsum(:,k)./(Ncount(k)+1);%( (mu0(k,:)./sigma0.^2) + ((Xsum(k,:))./sigma2) )* 1./( (1./sigma0.^2) + (D*(Ncount(k))./sigma2) );
        mu{iter}(k,:)=mvnrnd( m,s2 );
    end
end
end

