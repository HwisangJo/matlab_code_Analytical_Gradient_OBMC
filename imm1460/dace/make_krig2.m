clc;
close all
clear;

dm=2;
%% LHS
sample_n=10;
dm=2;
x=lhsdesign(sample_n,dm);
x1=-20+60*x(:,1); % -20~40
x2=15*x(:,2); % 0~15
x=[x1 x2];
% save('min_sample.mat','x1','x2','x')
save('krig10_sample.mat','x1','x2','x')
scatter(x(:,1),x(:,2))

% load('set_samp1');
load('min_sample');
% x1=x_B(:,1); % -20~40
% x2=x_B(:,2); % 0~15

%% Model
y1=(x2-(5.1/(4*pi^2))*x1.^2+5/pi*x1-6).^2+10*(1-1/(8*pi))*cos(x1)+10;
y2=5*cos(0.25*pi*x1).*x2+10*x1;
y3=5*sin(2*x1-0.5*pi).*x2-10*x2.^2;

x=[x1, x2];
y=[y1, y2, y3];

%% Kriging
theta0=1*ones(1,dm);
lb_theta=0.1*ones(1,dm);
ub_theta=10*ones(1,dm);

[dmodel1,perf1]=dacefit([x1,x2],y1,@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
[dmodel2,perf2]=dacefit([x1,x2],y2,@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
[dmodel3,perf3]=dacefit([x1,x2],y3,@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);

% save('dmodel_data_minsample.mat','dmodel1','dmodel2','dmodel3','x','y')
% save('krig10_model.mat','dmodel1','dmodel2','dmodel3','x','y')

cd('D:\Research\Research_JHS\MLE_krig')
save('krig10_model.mat','dmodel1','dmodel2','dmodel3','x','y')

%%


%%