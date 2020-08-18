clc
close all
clear
%% Sample
sample_n=13;
dm=2;

x=lhsdesign(sample_n,dm);

x1=-20+60*x(:,1); % -20~40 
x2=15*x(:,2); % 0~15

%% Model
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

X=gridsamp([-20 0;40 15],40);
result1=predictor(X,dmodel1);
result2=predictor(X,dmodel2);
result3=predictor(X,dmodel3);

X1=reshape(X(:,1),40,40);
X2=reshape(X(:,2),40,40);
result1=reshape(result1,40,40);
result2=reshape(result2,40,40);
result3=reshape(result3,40,40);
%% plot
figure(1)
mesh(X1,X2,result1)
figure(2)
contour(X1,X2,result1)
hold on
%%
p = -20:0.1:40;
q = 0:0.1:15;
figure(1)
scatter(x1,x2,1000,'k','.')
title('LHS 13 points')
xlabel('x1')
ylabel('x2')
hold on
scatter(Bay_mean1(1),Bay_mean1(2),1000,'k','+')
hold on
scatter(x_B(15,1),x_B(15,2),1000,'k','x')
hold on
scatter(x_B(16,1),x_B(16,2),1000,'k','x')
hold on
scatter(x_B(17,1),x_B(17,2),1000,'k','x')
legend('LHS','1st seq','2nd seq','Location','Best')
yy=CI(1,2)*ones(1,length(p));
hold on
plot(p,yy,'r','linewidth',1)
yy=CI(2,2)*ones(1,length(p));
hold on
plot(p,yy,'r','linewidth',1)
xx=CI(1,1)*ones(1,length(q));
hold on
plot(xx,q,'r','linewidth',1)
xx=CI(2,1)*ones(1,length(q));
hold on
plot(xx,q,'r','linewidth',1)
%%

figure(2)
scatter(x(3000:30000,1,2),x(3000:30000,2,2),'.')
hold on
scatter(x(3000:30000,1,3),x(3000:30000,2,3),'.')
hold on
scatter(x(3000:30000,1,4),x(3000:30000,2,4),'.')

figure(3)
mesh(X1,X2,result2)
figure(4)
contour(X1,X2,result2)
hold on
scatter(x1,x2,'k','.')

figure(5)
mesh(X1,X2,result3)
figure(6)
contour(X1,X2,result3)
hold on
scatter(x1,x2,'k','.')
%%