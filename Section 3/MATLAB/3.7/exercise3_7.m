clear all;close all;clc
rawData = readtable('tea_discipline_oss.csv','ReadVariableNames',true);
gradeCell = rawData.GRADE;
N = length(gradeCell);
action = rawData.ACTIONS;
grade = zeros(N,1);
delRow = [];
for index = 1:N
	try
		grade(index) = str2num(gradeCell{index});
	catch
		delRow = [delRow,index];
	end
end
action(delRow) = [];
grade(delRow) = [];
I = (action ~= -99);
action = action(I);
Y = action;
grade = grade(I);

X = [ones(length(grade),1),grade];
mu0 = [0;0];sigma0 = 1;
logPosBeta = @(beta) exp(-1/2/sigma0^2*(beta-mu0)'*(beta-mu0)) + sum(action.*(X*beta) - exp(X*beta));
negPos = @(beta) -logPosBeta(beta);
[betaPosMean,fval,exitflag,output] = fminsearch(negPos,[0.1;0.1]);
H11 = -1/sigma0^2 - sum(exp(X*betaPosMean) .*X(:,1).^2);
H22 = -1/sigma0^2 - sum(exp(X*betaPosMean) .*X(:,2).^2);
H12 = -1/sigma0^2 - sum(exp(X*betaPosMean) .*X(:,1).* X(:,2));
H = [H11,H12;H12,H22];
SIGMA = inv(-H);
CIbeta0 = [betaPosMean(1) - 1.96*sqrt(SIGMA(1,1)),betaPosMean(1) + 1.96*sqrt(SIGMA(1,1))]
CIbeta1 = [betaPosMean(2) - 1.96*sqrt(SIGMA(2,2)),betaPosMean(2) + 1.96*sqrt(SIGMA(2,2))]