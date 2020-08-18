clc;
close all
clear;

%% Sample
load('result_krig');

b_in=3000;

xx1=[];
xx2=[];
xx3=[];
xx4=[];

    j1=1;j2=1;j3=1;j4=1;
    for i=b_in:30000
        if imim(i,1)==1
            xx1(j1,1)=x(i,1,1);
            xx1(j1,2)=x(i,2,1);
            j1=j1+1;
        elseif imim(i,2)==1
            xx2(j2,1)=x(i,1,2);
            xx2(j2,2)=x(i,2,2);
            j2=j2+1;
        elseif imim(i,3)==1        
            xx3(j3,1)=x(i,1,3);
            xx3(j3,2)=x(i,2,3);
            j3=j3+1;
        elseif imim(i,4)==1            
            xx4(j4,1)=x(i,1,4);
            xx4(j4,2)=x(i,2,4);
            j4=j4+1;
        end
    end
%% Bayes 결과
xx=[xx2;xx3;xx4];
Bay_mean1=mean(xx);
Bay_var1=var(xx);
Bay_sd1=sqrt(Bay_var1);
% sum((xx(:,1)-Bay_mean1(1)).^2)/(length(xx)-1);
%% Bayes Confidence Interval
v=length(xx)-1;
td=-tinv(0.025,v)*1900;
CI=[Bay_mean1-td*Bay_sd1/sqrt(v+1);Bay_mean1+td*Bay_sd1/sqrt(v+1)];
%% 데이터 저장
%% New kriging model considering bayesian result
load('krig_data.mat')
load('Bayes_result.mat')

% x_B=[[x1 x2];Bay_mean1;x];
x_B=[[x1 x2];xopt];
x1=x_B(:,1);
x2=x_B(:,2);

% dm=2;
y1_B=(x2-(5.1/(4*pi^2))*x1.^2+5/pi*x1-6).^2+10*(1-1/(8*pi))*cos(x1)+10;
y2_B=5*cos(0.25*pi*x1).*x2+10*x1;
y3_B=5*sin(2*x1-0.5*pi).*x2-10*x2.^2;
y_B=[y1_B y2_B y3_B];

% theta0=1*ones(1,dm);
% lb_theta=0.1*ones(1,dm);
% ub_theta=10*ones(1,dm);
% 
% [dmodel1,perf1]=dacefit([x1,x2],y1,@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
% [dmodel2,perf2]=dacefit([x1,x2],y2,@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
% [dmodel3,perf3]=dacefit([x1,x2],y3,@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);

%% Kriging model of CVE
err=[];
for i=1:length(x_B)
x_B_=x_B;
y_B_=y_B;
x_B_(i,:)=[];
x1=x_B_(:,1);
x2=x_B_(:,2);

dm=2;
y_B_(i,:)=[];

theta0=1*ones(1,dm);
lb_theta=0.1*ones(1,dm);
ub_theta=10*ones(1,dm);

[dmodel1,perf1]=dacefit([x1,x2],y_B_(:,1),@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
[dmodel2,perf2]=dacefit([x1,x2],y_B_(:,2),@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
[dmodel3,perf3]=dacefit([x1,x2],y_B_(:,3),@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);

X=x_B(i,:);

result1=predictor(X,dmodel1);
result2=predictor(X,dmodel2);
result3=predictor(X,dmodel3);

y_i=[result1 result2 result3];
err(i,:)=y_B(i,:)-y_i;
end

x1=x_B(:,1);
x2=x_B(:,2);

dm=2;
theta0=1*ones(1,dm);
lb_theta=0.1*ones(1,dm);
ub_theta=10*ones(1,dm);

[dmodel1,perf1]=dacefit([x1,x2],err(:,1),@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
[dmodel2,perf2]=dacefit([x1,x2],err(:,2),@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);
[dmodel3,perf3]=dacefit([x1,x2],err(:,3),@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);

%% Optimization
x0=Bay_mean1-[1 1];
A = [];
b = [];
Aeq = [];
beq = [];

cve= @(x) -abs(min_dis(x)*cve_krig(x));
x_opt2=fmincon(cve,x0,A,b,Aeq,beq,CI(1,:),CI(2,:));


%% CVE graph
X=gridsamp([CI(1,:);CI(2,:)],40);
Y=[];
for i=1:length(X)
Y(i)=cve(X(i,:));
end
Y=Y';
X1=reshape(X(:,1),40,40);
X2=reshape(X(:,2),40,40);
YY=reshape(Y,40,40);


figure(1)
mesh(X1,X2,YY)
figure(2)
contour(X1,X2,YY)
%%


