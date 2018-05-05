%Morgan Kelley
%Ex. 5.2
clear all; close all; clc;
load restaurants.mat;

X=[restaurants.Profit];
X=Scale(X,0,1);
Y=restaurants.DinnerService;
Y=Scale(Y,0,1);

init=[1 1 .502 1 1 1];
Ns=100;


[mut,sigmat,omega1,zt,Pxi1,Pxi0,s2t,mt,Ncountt,Xcountt]=Ex52Func_2(Ns,Y,X,init);

error=sum(abs(Y'-zt(Ns-1,:)));