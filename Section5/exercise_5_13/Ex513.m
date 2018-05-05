%Morgan Kelley
%5.13
alpha=[0.01 0.1 1.0 10];
n=100;
ID=zeros(n,length(alpha));

ID1=cluster(alpha(1),n);
ID2=cluster(alpha(2),n);
ID3=cluster(alpha(3),n);
ID4=cluster(alpha(4),n);

SIZE=[length(ID1) length(ID2) length(ID3) length(ID4)];

cat=categorical({'\alpha_{0.01}','\alpha_{0.1}','\alpha_{1.0}','\alpha_{10}'});
bar(cat,SIZE)
ylabel('Number of clusters')