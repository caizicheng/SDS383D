function beta = gibbsSampler(X,Y,N)

	a0 = 0.1;
	b0 = 5;
	mu0 = zeros(size(X,2),1);
	K0 = 10 * diag(ones(1,size(X,2)));
	Kn = X'* X; + K0;
	an = a0 + 0.5 * size(X,1);
	A = ones(1,size(X,1));

	for index = 1:size(X,1)
		if Y(index) == 1
			A(index) = -1;
		end
	end

	w = ones(1,N);
	beta = ones(size(X,2),N);
	z = zeros(size(X,1),N);
	bn = ones(1,N) * b0;
	mu_n = ones(size(X,2),N);
	
	
	for index = 2:N
        index
		w(index) = gamrnd(an,bn(index - 1));
		bn(index) = b0; + 0.5*(z(:,index - 1)' * z(:,index - 1) - mu_n(:,index - 1)'*Kn*mu_n(:,index - 1));
		mu_n(:,index) = inv(Kn) * (X' * z(:,index - 1));
		try
			beta(:,index) = mvnrnd(mu_n(:,index),inv(w(index)  * Kn));
		catch
			keyboard
		end
		z(:,index) = rmvnrnd(X * beta(:,index),diag(ones(size(X,1),1) / sqrt(w(index))),1,A,0);
	
	end
	


end