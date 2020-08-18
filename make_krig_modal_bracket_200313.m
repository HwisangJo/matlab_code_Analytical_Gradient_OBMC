clc
close all
clear

%% Kriging model
sample_n=20;dm=2;

% Range
CI(1,:)=[norminv(0.0005,210e9,210e9*0.20) norminv(0.9995,210e9,210e9*0.20)];
CI(2,:)=[norminv(0.0005,8000,8000*0.20) norminv(0.9995,8000,8000*0.20)];

load('initial_samples_modal_bracket_200313.mat')

dm=2;
theta0=1*ones(1,dm);
lb_theta=0.1*ones(1,dm);
ub_theta=10*ones(1,dm);

addpath([pwd,'\imm1460\dace']) % Kriging toolbox
[dmodel1,perf1]=dacefit(x,y(:,1),@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);

dm=1;
theta0=1*ones(1,dm);
lb_theta=0.1*ones(1,dm);
ub_theta=10*ones(1,dm);
[dmodel2,perf2]=dacefit(x(:,1),y2(:,1),@regpoly2,@corrgauss,theta0,lb_theta,ub_theta);

save('dmodel_data_modal_bracket_200313.mat','dmodel1','dmodel2','x','y','y2') % 사전 분포 범위, 모르는 입력 분포 범위가 같을 때

%% Kriging Response - Graph
addpath([pwd,'\imm1460\dace']) % Kriging toolbox
CI(1,:)=[norminv(0.0005,210e9,210e9*0.20) norminv(0.9995,210e9,210e9*0.20)];
CI(2,:)=[norminv(0.0005,8000,8000*0.20) norminv(0.9995,8000,8000*0.20)];

Fig=figure;
set(Fig,'pos',[700 541 710 553]);
[X,Y] = meshgrid(linspace(CI(1,1),CI(1,2),100),linspace(CI(2,1),CI(2,2),100));
for i=1:size(X,1)
    Z(:,i)=model_krig_modal_bracket([X(:,i) Y(:,i)],2);
end
surf(X,Y,Z,'FaceAlpha',0.5,'EdgeColor','none')
xlim(CI(1,:))
ylim(CI(2,:))
box on
hold on
scatter3(x(:,1),x(:,2),y,50,'MarkerEdgeColor','w','MarkerFaceColor',[1 0 0])
xlabel('Youngs Modulus (\theta_{1})')
ylabel('Density (\theta_{2})')
set(get(gca,'XLabel'),'Rotation',15)
set(get(gca,'YLabel'),'Rotation',-25)
zlabel('1st Natural Frequency (y_{2})')
grid on
legend({'Kriging Model','Simulations'},'location','best')
set(gca,'fontname','times')  % Set it to times
axx = gca;
axx.FontSize = 15;

figure
x_i=linspace(CI(1,1),CI(1,2),50)';
% x_i=linspace(0,CI(2,2),50)';
y_i=-model_krig_modal_bracket2(x_i,1);
plot(x_i,y_i,'k-','linewidth',2.0)
xlim(CI(1,:))
box on
hold on
scatter(x(:,1),-y2(:,1),50,'MarkerEdgeColor','w','MarkerFaceColor',[1 0 0])
xlabel('Youngs Modulus (\theta_{1})')
ylabel('Deflection (y_{1})')
grid on
legend({'Kriging Model','Simulations'})
set(gca,'fontname','times')  % Set it to times
axx = gca;
axx.FontSize = 15;
