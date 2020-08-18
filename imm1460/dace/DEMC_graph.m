
clc;clear;

%% DE-MC Marginal distribution

plot(1:t,x(:,1,1))
hist(x(:,1,1),50)
hist(x(:,1,2),50)
hist(x(:,1,3),50)
hist(x(:,1,4),50)

figure(2)
k=[x(:,1,1);x(:,1,2);x(:,1,3);x(:,1,4)];
hist(k,100)
xlabel('ks')
ylabel('N')
title('Marginal pdf hist(DE-MC Sampling)')

sim2=k;
pdsim2 = fitdist(sim2,'Kernel')
x4 = 10000:50:25000;
yyy = pdf(pdsim2,x4);
yyyy=yyy/(sum(yyy)*50);
figure(10)
plot(x4,yyyy)
xlabel('ks')
ylabel('pdf of ks')
title('Marginal pdf(DE-MC Sampling)')

figure(3)
k=[x(:,2,1);x(:,2,2);x(:,2,3);x(:,2,4)];
hist(k,100)
xlabel('kt')
ylabel('N')
title('Marginal pdf hist(DE-MC Sampling)')

figure(3)
plot(1:t,x(:,1,1))
xlabel('T')
ylabel('ks')
title('Search(DE-MC Sampling)')

%% DE-MC 2-variate sampling
figure(1)
scatter(xx(:,1),x(:,2),'.')
hold on
scatter(x(3000:30000,1,1),x(3000:30000,2,1),'.')
hold on
scatter(x(3000:30000,1,2),x(3000:30000,2,2),'.')
hold on
scatter(x(3000:30000,1,3),x(3000:30000,2,3),'.')
hold on
scatter(x(3000:30000,1,4),x(3000:30000,2,4),'.')
xlabel('ks')
ylabel('kt')
title('Joint pdf(DE-MC Sampling)')
%% DE-MC kernel density estimation
d=2;
n=length(x(:,1,1));

%covariance
covF=cov(x(:,1,1), x(:,2,1));
costdF=sqrt(covF);

%contour
p = 500:20:1500;
q = 15960:1.6:16040;
[P,Q] = meshgrid(p,q);

%bandwidth
HrtF = (4/(d+2))^(1/(d+4))*(n)^(-1/(d+4))*covF;

%Multivariate KDE
su=zeros(length(p),length(q));
for i=1:n
F = mvnpdf([P(:) Q(:)],[x(i,1,1) x(i,2,1)],HrtF);
F = reshape(F,length(q),length(p));
su(:,:) = su(:,:) + F(:,:);
end
su = su/n;

figure(9)
contour(P,Q,su)
ylabel('ks');
xlabel('bs');
title('Joint pdf of ks and bs(DE-MC)');
%colormap('jet')
colorbar

figure(10)
plot(x_check(1:1000,1,1),x_check(1:1000,2,1),'+')
ylabel('ks');
xlabel('bs');
title('sampling of bs and ks(DE-MC)');

hold on
plot(x(100:1000,1,1),x(100:1000,2,1),'o')
ylabel('ks');
xlabel('bs');
title('sampling of bs and ks(DE-MC)');

% 참값, 수직선 plot
yy=16000*ones(1,length(q));
hold on
plot(p,yy,'k','linewidth',1)

xx=300*ones(1,length(p));
hold on
plot(xx,q,'k','linewidth',1)

% 보정된 값
calib_x=mean(x(18000:30000,1,1));
calib_y=mean(x(18000:30000,2,1));

yy=calib_y*ones(1,length(q));
hold on
plot(p,yy,'r','linewidth',1)

xx=calib_x*ones(1,length(p));
hold on
plot(xx,q,'r','linewidth',1)

% 편차
dev_x=calib_x-300;
dev_y=calib_y-16000;

hold on
plot(x(500:30000,1,1),x(500:30000,2,1),'+')
ylabel('ks');
xlabel('mb');
title('sampling of mb and ks(DE-MC)');

%% DE-MC 움직임
plot(x(1:20,1,1),x(1:20,2,1))
hold on
plot(x(1:20,1,2),x(1:20,2,2))
hold on
plot(x(1:20,1,3),x(1:20,2,3))
hold on
plot(x(1:20,1,4),x(1:20,2,4))
ylabel('ks');
xlabel('mw');
title('DE-MC sampling itr=20');
%% RWM
clc
close all
clear

origin=[12800 190000];
A=origin*0.3;
B=2*origin;

pdf= @(x) object9_var2(x);
prior = @(N,d) unifrnd(A,B);

T=30000;
d=2;

q = @(C,d) mvnrnd(zeros(1,d),C); % Multivariate normal proposal distribution
C = (2.38/sqrt(d))^2 * eye(d); % Covariance matrix proposal distribution
x = nan(T,d); p_x = nan(T,1); % Preallocate memory for chain and density
x(1,1:d) = prior(1,d); % Initialize chain by sampling from prior
p_x(1) = pdf(x(1,1:d)); % Calculate density current state chain

for t = 2:T, % Dynamic part: Chain evolution
xp = x(t-1,1:d) + q(C,d); % Create candidate point
p_xp = pdf(xp); % Calculate density of proposal
alpha = min(p_xp/p_x(t-1),1); % Compute Metropolis ratio
idx = alpha > rand; % Alpha larger than U[0,1] or not?
if idx, % Idx = 0 (false) or 1 (true)
x(t,1:d) = xp; p_x(t) = p_xp; % True: accept proposal
else
x(t,1:d) = x(t-1,1:d); p_x(t) = p_x(t-1); % False: stay at current state of chain
end
end

% for t = 2:T, % Dynamic part: Chain evolution
% xp = x(t-1,1:d) + q(C,d); % Create candidate point
% p_xp = pdf(xp); % Calculate density of proposal
% ttt=p_xp-p_x(t-1);
% alpha = min(ttt,0); % Compute Metropolis ratio
% idx = alpha > log(rand); % Alpha larger than U[0,1] or not?
% if idx, % Idx = 0 (false) or 1 (true)
% x(t,1:d) = xp; p_x(t) = p_xp; % True: accept proposal
% else
% x(t,1:d) = x(t-1,1:d); p_x(t) = p_x(t-1); % False: stay at current state of chain
% end
% end
%% RWM graph

sim2=x(10000:30000,1,1);
pdsim2 = fitdist(sim2,'Kernel');
x4 = 10000:50:30000;
yyy = pdf(pdsim2,x4);

figure
plot(x4,400*yyy)
figure
histogram(sim2, 'Normalization' , 'probability')

scatter(x(:,1),p_x,'.')
scatter(x(:,2),p_x,'.')
scatter(x(:,1),x(:,2),'.')
hist(x(:,1),50)

%%