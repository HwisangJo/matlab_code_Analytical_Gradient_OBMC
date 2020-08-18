clc
close all
clear

%% Input
% x1=normrnd(2,0.1^2,[5,1]);
x1=wblrnd(5.7,1.3,[10,1]);
% x2=normrnd(8,0.2^2,[5,1]);
x2=lognrnd(2,0.3,[10,1]);
x=[x1 x2];
%% Model
y1=(x2-(5.1/(4*pi^2))*x1.^2+5/pi*x1-6).^2+10*(1-1/(8*pi))*cos(x1)+10;
y2=5*cos(0.25*pi*x1).*x2+10*x1;
y3=5*sin(2*x1-0.5*pi).*x2-10*x2.^2;

exp_data=[y1 y2 y3];
%% Sample mena, variance
mean_samp=mean(x);
var_samp=var(x);
sd_samp=std(x);

%% ¿˙¿Â
save('exp_data2.mat','exp_data','mean_samp','sd_samp','var_samp','x','x1','x2','y1','y2','y3')

%%