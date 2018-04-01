function [beta] = Bcond_func(tau,Bn)
%Conditional beta generator
% Kn[8x8], Bn[1x8], tau[1x1]
global Kn;
sigma=(Kn*tau)^-1;
mu=Bn;
beta=mvnrnd(mu,sigma);
end

