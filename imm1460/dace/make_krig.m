clc;
close all
clear;

sample_n=37;
% sample_n=24;
dm=2;
%% LHS
x=lhsdesign(sample_n,dm);
x1=-20+60*x(:,1); % -20~40
x2=15*x(:,2); % 0~15
x=[x1 x2];
save('37_sample.mat','x1','x2','x','dmodel1','dmodel2','dmodel3')
save('min_sample.mat','x1','x2','x')
scatter(x(:,1),x(:,2))
% xa=lhsdesign(sample_n,dm);
% xa1=-20+60*xa(:,1); % -20~40
% xa2=15*xa(:,2); % 0~15
% 
% x=[x;[xa1 xa2]];
%% Sample ±×·¡ÇÁ Ç¥½Ã
clc
close all
clear

load('krig_data_original.mat')
scatter(x(:,1),x(:,2),90,'k','filled')
p = 1:0.1:3.5;
q = 5:0.1:10;
yy=5*ones(1,length(p));
hold on
plot(p,yy,'k','linewidth',2)
xx=1*ones(1,length(q));
hold on
plot(xx,q,'k','linewidth',2)
yy=10*ones(1,length(p));
hold on
plot(p,yy,'k','linewidth',2)
xx=3.5*ones(1,length(q));
hold on
plot(xx,q,'k','linewidth',2)
ylabel('x2');
xlabel('x1');
title('LHS 13 points');

% Å©¸®±ë »ùÇÃ & ½ÇÇè°ª »ý¼º ÀÎÇ² »ùÇÃ Plot
clc
close all
clear

load('min_sample.mat');
scatter(x(:,1),x(:,2),100,'b')
ylabel('x2');
xlabel('x1');
title('1st Kriging Model');
load('exp_data.mat');
hold on
scatter(x(:,1),x(:,2),30,'k','filled')
legend('Kriging Sample','Experiment Input Sample','Location','Best');

% Å©¸®±ë »ùÇÃ & Sequential Sample Plot
load('min_sample');
% load('result_seq_samp1.mat')x(:,1),x(:,2)
scatter(x(:,1),x(:,2),90,'k','filled')
% scatter(x(1:13,1),x(1:13,2),90,'k','filled')
hold on
scatter(x(14:19,1),x(14:19,2),350,'r','x')
hold on
scatter(1.9993,8.0312,300,'b','+')
% legend('LHS','Seq.Sample1','Actual Input','Location','Best');
ylabel('x2');
xlabel('x1');
title('Sequential Sampling');

p = CI(1,1):0.001:CI(1,2);
q = CI(2,1):0.001:CI(2,2);
yy=CI(2,1)*ones(1,length(p));
hold on
plot(p,yy,'r','linewidth',2)
xx=CI(1,1)*ones(1,length(q));
hold on
plot(xx,q,'r','linewidth',2)

yy=CI(2,2)*ones(1,length(p));
hold on
plot(p,yy,'r','linewidth',2)
xx=CI(1,2)*ones(1,length(q));
hold on
plot(xx,q,'r','linewidth',2)
axis([1 3.5 5 10])

load('result_seq_samp2.mat')
hold on
scatter(x(20:22,1),x(20:22,2),350,'g','x')
p = CI(1,1):0.001:CI(1,2);
q = CI(2,1):0.001:CI(2,2);
yy=CI(2,1)*ones(1,length(p));
hold on
plot(p,yy,'g','linewidth',2)
xx=CI(1,1)*ones(1,length(q));
hold on
plot(xx,q,'g','linewidth',2)

yy=CI(2,2)*ones(1,length(p));
hold on
plot(p,yy,'g','linewidth',2)
xx=CI(1,2)*ones(1,length(q));
hold on
plot(xx,q,'g','linewidth',2)
axis([1 3.5 5 10])

load('result_seq_samp3.mat')
hold on
scatter(x(23,1),x(23,2),350,'m','x')
p = CI(1,1):0.001:CI(1,2);
q = CI(2,1):0.001:CI(2,2);
yy=CI(2,1)*ones(1,length(p));
hold on
plot(p,yy,'m','linewidth',2)
xx=CI(1,1)*ones(1,length(q));
hold on
plot(xx,q,'m','linewidth',2)

yy=CI(2,2)*ones(1,length(p));
hold on
plot(p,yy,'m','linewidth',2)
xx=CI(1,2)*ones(1,length(q));
hold on
plot(xx,q,'m','linewidth',2)
axis([1 3.5 5 10])


%% sample È®ÀÎ
scatter(x(:,1),x(:,2),'k','.')
p = -20:0.1:40;
q = 0:0.1:15;
yy=8*ones(1,length(p));
hold on
plot(p,yy,'r','linewidth',1)
xx=2*ones(1,length(q));
hold on
plot(xx,q,'r','linewidth',1)

%% Model
x1=x(:,1);x2=x(:,2);
y1=(x2-(5.1/(4*pi^2))*x1.^2+5/pi*x1-6).^2+10*(1-1/(8*pi))*cos(x1)+10;
y2=5*cos(0.25*pi*x1).*x2+10*x1;
y3=5*sin(2*x1-0.5*pi).*x2-10*x2.^2;

%% Kriging
theta0=1*ones(1,dm);
lb_theta=0.1*ones(1,dm);
ub_theta=10*ones(1,dm);

[dmodel1,perf1]=dacefit([x1,x2],y1,@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
[dmodel2,perf2]=dacefit([x1,x2],y2,@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
[dmodel3,perf3]=dacefit([x1,x2],y3,@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
%%