function J = objFun(y,f,K)
	
	p_y_f = @(yi,fi) (1/(1+exp(-fi)))^yi * (1-1/(1+exp(-fi)))^(1-yi);
	sumLogPfy = 0;
	for index = 1:length(y)
		sumLogPfy = sumLogPfy + log(p_y_f(y(index),f(index)));
	end
	sumLogPfy = sumLogPfy - 0.5 * f'*inv(K)*f - 0.5*log(det(K));
	J = -sumLogPfy;
end