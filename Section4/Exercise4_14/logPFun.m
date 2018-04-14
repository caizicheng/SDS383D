function logP = logPFun(PhiX,y,alpha2,l2,sigma2,x)
	
	
	K_gen = @(x1,x2) alpha2*exp(-1/2/l2 * norm(x1-x2,2)^2);
	
	K__ = zeros(length(x),length(x));

	for row = 1:length(x)
		for col = 1:length(x)
			K__(row,col) = K_gen(PhiX(row,:),PhiX(col,:)) + sigma2*(row == col);
		end
	end
	
	logP = -0.5*y'*inv(K__)*y - 0.5*log(det(K__));

end