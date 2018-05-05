function [mu,z,Ncount] = Ex57Func_7(Ns,X,K,mu0,sigma0)
[N,D] = size(X);
n=N/K;
%Initial values

%omega=zeros(Ns+1,K);
omega = 1;
%Initialize parameter size
mu=cell(Ns+1);%;zeros(Ns+1,d);
z=ones(Ns+1,N);
z(1,:)=[5*ones(1,n) 3*ones(1,n) 4*ones(1,n) ones(1,n) 2*ones(1,n)];

alpha=1;
s2=eye(D);
sigma2=sigma0;
%Initial Mu0
mu{1} = mvnrnd(mu0,eye(D));

Ncount=zeros(1,K);
Ncount(1) = N;
Xsum = zeros(K,D);
Xsum(1,:) = sum(X,1);
r=rand(N,1);

for iter=2:Ns+1
   a=1;
    for i=1:N
        
        for k=1:K
            Mcount(i,k)=Ncount(k) -double(z(iter-1,i)==k)+ alpha;
            
            prob(i,k) = Mcount(i,k)*mvnpdf(X(i,:),mu{iter-1}(k,:),eye(D));
        end
        
        prob(i,:) = prob(i,:)/sum(prob(i,:));

        z(iter,i)= sample2(prob(i,:),1);

   
    end
    
    Xsum = [sum(X(z(iter,:)==1,:),1)' sum(X(z(iter,:)==2,:),1)' sum(X(z(iter,:)==3,:),1)' sum(X(z(iter,:)==4,:),1)' sum(X(z(iter,:)==5,:),1)'];
    
    Ncount = [sum(z(iter,:)==1) sum(z(iter,:)==2) sum(z(iter,:)==3) sum(z(iter,:)==4) sum(z(iter,:)==5)];   
    sigma2=1;%/omega;
    mu{iter} = zeros(K,D);
    
    for k=1:K
        s2=eye(2)./(Ncount(k)+1);%( (1./sigma0.^2) + (D*(Ncount(k))./sigma2) ).^(-1);
        m=Xsum(:,k)./(Ncount(k)+1);%( (mu0(k,:)./sigma0.^2) + ((Xsum(k,:))./sigma2) )* 1./( (1./sigma0.^2) + (D*(Ncount(k))./sigma2) );
        mu{iter}(k,:)=mvnrnd( m,s2 );
    end
    
end
end

