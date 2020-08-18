%% MLE inverse problem : KD-estimation
clc
close all
clear

folder_path=pwd;
cd(folder_path)

opt_folder=[pwd,'\opt_data_ex2_krig_200309'];
mkdir(opt_folder)

addpath(genpath([pwd,'\imm1460\dace']))
% folder_name_type={'\exp_data_case1_normal_gumbel_krig_200309';
%     '\exp_data_case2_lognormal_extreme_krig_200309';
%     };
% dist_type={'Normal','Extreme Value';'Lognormal','Extreme Value'};
% dist_name_type={'case1_normal_gumbel';'case2_lognormal_extreme'};
% type_num_case=[1 3;2 6];
folder_name_type={'\exp_data_case1_normal_lognormal_krig_200309';
    '\exp_data_case2_gumbel_extreme_krig_200309';
    };
dist_type={'Normal','Lognormal';'Extreme Value','Extreme Value'};
dist_name_type={'case1_normal_lognormal';'case2_gumbel_extreme'};
type_num_case=[1 2;3 6];
% dpt_n=[10 8]; % dimension and number of experiment results
dpt_n=[50 1]; % dimension and number of experiment results
% dpt_n=[200 1]; % dimension and number of experiment results
%% Test - Input variable
% 시험 데이터 생성, disk, 190414
for nt=1:2
    % 저장
        folder_name='\exp_data_ex2_bracket_200316';
file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
    sprintf(['exp_data_ex2_krig_',dist_name_type{nt},'_n%d_%d_trial%d_200309.mat'],dpt_n(2),dpt_n(1),1)];
    load(file_exp_name)
    for ii=1:2
        x_sample(:,ii,nt)=reshape(xc(:,ii,:),size(xc,1)*size(xc,3),1);
    end
end

% type_num_case=[1 3;2 6];
% dist_type={'Normal','Extreme Value';'Lognormal','Extreme Value'};
for nt=1:2
    for ii=1:2
    phat1 = mle(x_sample(:,ii,nt),'distribution',dist_type{nt,ii});
    phat2=phat1;
    if type_num_case(nt,ii)==2
        [M,N]=lognstat(phat1(1),phat1(2));
        phat2=[M,sqrt(N)];
    elseif type_num_case(nt,ii)==3
        phat1 = mle(-x_sample(:,ii,nt),'distribution',dist_type{nt,ii});
%         [M,N]=evstat(phat1(1),phat1(2));
        phat2=[-phat1(1)+0.5772*phat1(2),sqrt((-phat1(2))^2*pi^2/6)];
    elseif type_num_case(nt,ii)==6
        [M,N]=evstat(phat1(1),phat1(2));
        phat2=[M,sqrt(N)];
    elseif type_num_case(nt,ii)==5
        [M,N]=unifstat(phat1(1),phat1(2));
        phat2=[M,sqrt(N)];
    end
    parm_true2(nt,[2*(ii-1)+1:2*(ii-1)+2])=phat2;
    end
end
% parm_true2(:,:,1);
foldername4= [opt_folder,'\','opt_input_mle'];
save(foldername4,'parm_true2') % sim?
%% Contour Plot
foldername5=[opt_folder,'\','opt_hist_analy'];
mkdir(foldername5)
for nt=2%:2
load([pwd,'\opt_data_ex2_krig_200309\opt_input_mle.mat'])
    folder_name='\exp_data_ex2_bracket_200316';
    for ini=2%:2%:5
for i=1%:6
foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_norm_trr_analy.mat', nt,ini,i)];

% Optimization-Based Statistical Model Calibration (Parameter Estimation)
file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
    sprintf(['exp_data_ex2_krig_',dist_name_type{nt},'_n%d_%d_trial%d_200309.mat'],dpt_n(2),dpt_n(1),1)];
file_s_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',i)];
save_sign=0;

% Performance Function
perf_func1=@(x) model_krig_modal_bracket(x,2);
perf_func2=@(x) model_krig_modal_bracket2(x,1);
perf_func=@(x,y_num) model_krig_mult_resp(x,y_num,perf_func1,perf_func2);

lbx=[210e9, 210e9*0.10, 8000, 8000*0.10]*0.5;
ubx=[210e9, 210e9*0.10, 8000, 8000*0.10]*1.5;

ratio_x0=[1.10 1.05 0.90 0.95;0.92 0.98 1.10 1.10; 0.95 0.95 1.10 1.05;1.03 1.07 0.98 0.92;1.04 0.95 0.92 0.90];
x0=parm_true2(nt,:).*ratio_x0(ini,:);

load(foldername3)
end

    end
end

% scatter(y_exp([1:150],1),y_exp([151:300],1),'k.')
% scatter(v_exp_data(1,1,1:50),v_exp_data(1,2,1:50),'k.')
CI(1,:)=[norminv(0.0005,210e9,210e9*0.10) norminv(0.9995,210e9,210e9*0.10)];
CI(2,:)=[norminv(0.0005,8000,8000*0.10) norminv(0.9995,8000,8000*0.10)];

% y_range1=[120 200];
% y_range2=[-10 -7.0]*10^-3;
% y_range1_=[120 200];
% y_range2_=[-10 -7.0]*10^-3;

% y_range1=[150 170];
% y_range2=[-7.5 -9.5]*10^-3;
% y_range1_=[165 190];
% y_range2_=[-8.5 -7.0]*10^-3;

% y_range1_=[120 165];
% y_range2_=[-7.5 -12.0]*10^-3;
% y_range1=[140 180];
% y_range2=[-7.0 -10.5]*10^-3;

y_range1_=[120 165];
y_range2_=[7.5 12.0]*10^-3;
y_range1=[140 180];
y_range2=[7.0 10.5]*10^-3;

load(file_s_name)
load(file_exp_name)
n=10^6;
sim=zeros(n,1);
rng(s);
dim=size(x0,2);
X_=x0;
% type_num_=type_num_case(1,:);
type_num_=type_num_case(2,:);
dist_type_={'Normal','Lognormal','Extreme Value',[],'Uniform','Extreme Value'};
dim=size(X_,2);
for h=1:dim/2
    if type_num_(h)==5
        ab(1:2)=convert_param(X_((h-1)*2+1),X_(h*2),5);
        X((h-1)*2+1)=ab(1)-ab(2)*sqrt(3);
        X(h*2)=ab(1)+ab(2)*sqrt(3);
    else
        X((h-1)*2+1:h*2)=convert_param(X_((h-1)*2+1),X_(h*2),type_num_(h));
    end
end
for h=1:dim/2
    if type_num_(:,h)==3
        var_all(:,h)=X((h-1)*2+1)-(random(dist_type_{type_num_(:,h)},X((h-1)*2+1),X(h*2), [n, 1])-X((h-1)*2+1)); % Unknown Variables
    else
        var_all(:,h)=random(dist_type_{type_num_(:,h)},X((h-1)*2+1),X(h*2), [n, 1]); % Unknown Variables
    end
end
for i=1:1:size(v_exp_data,2) %y
%     Xd=xd(:,i);
    Xd=[];
%     parfor j=1%:size(xd,1)
    for j=1%:size(xd,1)
%         y=perf_func([repmat(Xd(j,:),n,1), var_all],i);
        y=perf_func([var_all],i);
        simm(:,j,i)=y;
        if i==2
        simm(:,j,i)=-y;
        end
    end
end

d=size(v_exp_data,2);
m=(4/((d+2)*n))^(1/(d+4));
j=1;
for i=1:1:size(v_exp_data,2) %y
    mean_y(j,i)=mean(simm(:,j,i));
    sigma_y(j,i)=sqrt(sum((simm(:,j,i)-mean_y(j,i)).^2/(n-1)));
    %         bw(j,i)=(4/(3*n))^(1/5)*sigma_y(j,i);
    bw(j,i)=m*sigma_y(j,i);
end

[X1,Y1] = meshgrid(linspace(y_range1_(1),y_range1_(2),50),linspace(y_range2_(1),y_range2_(2),50));
for k=1:size(X1,1)
    grid_=[X1(:,k) Y1(:,k)];
    for t=1:size(X1,1)
        for i=1:1:size(v_exp_data,2) %y
            pd_1(:,t,i)=normpdf(grid_(t,i),simm(:,1,i),bw(1,i));
        end
    end
    pd__(k,:)=sum(prod(pd_1,3),1)/size(simm,1);
    disp(sprintf('%d번째 완료',k))
end

n=10^6;
sim=zeros(n,1);
rng(s);
dim=size(x0,2);
X_=xvar;
% type_num_=type_num_case(1,:);
type_num_=type_num_case(2,:);
dist_type_={'Normal','Lognormal','Extreme Value',[],'Uniform','Extreme Value'};
dim=size(X_,2);
for h=1:dim/2
    if type_num_(h)==5
        ab(1:2)=convert_param(X_((h-1)*2+1),X_(h*2),5);
        X((h-1)*2+1)=ab(1)-ab(2)*sqrt(3);
        X(h*2)=ab(1)+ab(2)*sqrt(3);
    else
        X((h-1)*2+1:h*2)=convert_param(X_((h-1)*2+1),X_(h*2),type_num_(h));
    end
end
for h=1:dim/2
    if type_num_(:,h)==3
        var_all(:,h)=X((h-1)*2+1)-(random(dist_type_{type_num_(:,h)},X((h-1)*2+1),X(h*2), [n, 1])-X((h-1)*2+1)); % Unknown Variables
    else
        var_all(:,h)=random(dist_type_{type_num_(:,h)},X((h-1)*2+1),X(h*2), [n, 1]); % Unknown Variables
    end
end
for i=1:1:size(v_exp_data,2) %y
%     Xd=xd(:,i);
    Xd=[];
%     parfor j=1%:size(xd,1)
    for j=1%:size(xd,1)
%         y=perf_func([repmat(Xd(j,:),n,1), var_all],i);
        y=perf_func([var_all],i);
        simm(:,j,i)=y;
        if i==2
        simm(:,j,i)=-y;
        end
    end
end

d=size(v_exp_data,2);
m=(4/((d+2)*n))^(1/(d+4));
j=1;
for i=1:1:size(v_exp_data,2) %y
    mean_y(j,i)=mean(simm(:,j,i));
    sigma_y(j,i)=sqrt(sum((simm(:,j,i)-mean_y(j,i)).^2/(n-1)));
    %         bw(j,i)=(4/(3*n))^(1/5)*sigma_y(j,i);
    bw(j,i)=m*sigma_y(j,i);
end

[X2,Y2] = meshgrid(linspace(y_range1(1),y_range1(2),50),linspace(y_range2(1),y_range2(2),50));
for k=1:size(X2,1)
    grid_=[X2(:,k) Y2(:,k)];
    for t=1:size(X2,1)
        for i=1:1:size(v_exp_data,2) %y
            pd_1(:,t,i)=normpdf(grid_(t,i),simm(:,1,i),bw(1,i));
        end
    end
    pd_(k,:)=sum(prod(pd_1,3),1)/size(simm,1);
    disp(sprintf('%d번째 완료',k))
end

% [f,xi]=ksdensity([simm(:,1,1) simm(:,1,2)],[X(:),Y(:)],'Bandwidth',bw);
% ff=reshape(f,50,50);

Fig=figure;
set(Fig,'pos',[700 541 710 593]);
hold on
s2=contour(Y1,X1,pd__',':','LineWidth',2.0);
hold on
s1=contour(Y2,X2,pd_','-','LineWidth',1.5);
% hold on
% s2=contour(X,Y,ff,':','LineWidth',2.0);
hold on
% scatter(v_exp_data(1,1,1:50),v_exp_data(1,2,1:50),'rx')
% s3=scatter(v_exp_data(1,1,1:50),-v_exp_data(1,2,1:50),100,'k.');
s3=scatter(-v_exp_data(1,2,1:50),v_exp_data(1,1,1:50),100,'k.');
box on
% xlim([min(simm(:,1,1)),max(simm(:,1,1))])
% ylim([min(simm(:,1,2)),max(simm(:,1,2))])
% ylim([120 210])
% xlim([6.5*10^-3 11.5*10^-3])
ylim([120 210])
xlim([6.0*10^-3 11.5*10^-3])
axx = gca;
axx.FontSize = 17;
grid on
xlabel('Deflection (m) - \it y_{1}','fontsize',20)
ylabel('1st natural frequency (Hz) - \it y_{2}','fontsize',20)
grid on
legend({'Initial output PDF','Calibrated output PDF','Experimental data'},'Location','best','fontsize',17)
set(gca,'fontname','times')  % Set it to times


%% Estimated Input 1D
for nt=1:2
    % 저장
        folder_name='\exp_data_ex2_bracket_200316';
file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
    sprintf(['exp_data_ex2_krig_',dist_name_type{nt},'_n%d_%d_trial%d_200309.mat'],dpt_n(2),dpt_n(1),1)];
    load(file_exp_name)
    for ii=1:2
        x_sample(:,ii,nt)=reshape(xc(:,ii,:),size(xc,1)*size(xc,3),1);
    end
end

for nt=2%1:2%:2
    for ini=2
foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_norm_trr_analy.mat', nt,ini,1)];
load(foldername3)

type_num_=type_num_case(nt,:);
dist_type_={'Normal','Lognormal','Extreme Value',[],'Uniform','Extreme Value'};
X_=xvar;
dim=size(X_,2);
for h=1:dim/2
    if type_num_(h)==5
        ab(1:2)=convert_param(X_((h-1)*2+1),X_(h*2),5);
        X((h-1)*2+1)=ab(1)-ab(2)*sqrt(3);
        X(h*2)=ab(1)+ab(2)*sqrt(3);
    else
        X((h-1)*2+1:h*2)=convert_param(X_((h-1)*2+1),X_(h*2),type_num_(h));
    end
end

ratio_x0=[1.10 1.05 0.90 0.95;0.92 0.98 1.10 1.10; 0.95 0.95 1.10 1.05;1.03 1.07 0.98 0.92;1.04 0.95 0.92 0.90];
X0=parm_true2(nt,:).*ratio_x0(ini,:);
for h=1:dim/2
    if type_num_(h)==5
        ab(1:2)=convert_param(X0((h-1)*2+1),X0(h*2),5);
        X0_((h-1)*2+1)=ab(1)-ab(2)*sqrt(3);
        X0_(h*2)=ab(1)+ab(2)*sqrt(3);
    else
        X0_((h-1)*2+1:h*2)=convert_param(X0((h-1)*2+1),X0(h*2),type_num_(h));
    end
end

CI(1,1)=norminv(0.000005,xvar(1),xvar(2));CI(1,2)=norminv(0.999995,xvar(1),xvar(2));
CI(2,1)=norminv(0.000005,xvar(3),xvar(4));CI(2,2)=norminv(0.999995,xvar(3),xvar(4));
CI_(1,1)=norminv(0.000005,X0(1),X0(2));CI_(1,2)=norminv(0.999995,X0(1),X0(2));
CI_(2,1)=norminv(0.000005,X0(3),X0(4));CI_(2,2)=norminv(0.999995,X0(3),X0(4));
sample_n=100;
% x1 = linspace(min(CI(1,1),CI_(1,1)), max(CI(1,2),CI_(1,2)),sample_n);
% x2 = linspace(min(CI(2,1),CI_(2,1)), max(CI(2,2),CI_(2,2)),sample_n);
x1 = linspace(min(CI(1,1),CI_(1,1))+0.2*10^11, max(CI(1,2),CI_(1,2))+0.2*10^11,sample_n);
x2 = linspace(min(CI(2,1),CI_(2,1))-1000, max(CI(2,2),CI_(2,2))-1000,sample_n);
xx=[x1' x2'];

for iii=1:2
ypdf_ini(:,iii)=pdf(dist_type_{type_num_case(nt,iii)},xx(:,iii),X0_((iii-1)*2+1),X0_((iii-1)*2+2));
ypdf_cal(:,iii)=pdf(dist_type_{type_num_case(nt,iii)},xx(:,iii),X((iii-1)*2+1),X((iii-1)*2+2));
ypdf_mle(:,iii)=pdf(dist_type_{type_num_case(nt,iii)},xx(:,iii),parm_true2(nt,(iii-1)*2+1),parm_true2(nt,(iii-1)*2+2));
if type_num_case(nt,iii)==3
    ypdf_ini(:,iii)=pdf(dist_type_{type_num_case(nt,iii)},-xx(:,iii)+2*X0_((iii-1)*2+1),X0_((iii-1)*2+1),X0_((iii-1)*2+2));
    ypdf_cal(:,iii)=pdf(dist_type_{type_num_case(nt,iii)},-xx(:,iii)+2*X((iii-1)*2+1),X((iii-1)*2+1),X((iii-1)*2+2));
    ypdf_mle(:,iii)=pdf(dist_type_{type_num_case(nt,iii)},-xx(:,iii)+2*parm_true2(nt,(iii-1)*2+1),parm_true2(nt,(iii-1)*2+1),parm_true2(nt,(iii-1)*2+2));
end
end

nbin1=8;
Fig=figure(2*nt-1+10);
pp1=histogram(x_sample(:,1,nt),nbin1,'Normalization','pdf');
hold on
pp2=plot(x1,ypdf_ini(:,1),'r-.','LineWidth',2.0);
hold on
pp3=plot(x1,ypdf_cal(:,1),'k-','LineWidth',2.0);
xlim([min(CI(1,1),CI_(1,1))+0.2*10^11, max(CI(1,2),CI_(1,2))+0.2*10^11])
% ylim([0 6]*10^-5)
set(gca,'fontsize', 13,'fontname','times');
ytickangle(90)
xlabel('Youngs Modulus (Pa) - \theta_{1}','fontsize', 15)
ylabel('PDF','fontsize', 15)
legend([pp1 pp2 pp3],'Virtual input data','Prior PDF','Estimated PDF','Location','best','fontsize', 12)
set(Fig,'pos',[100 541 660 303]);

nbin2=8;
% nbin2=nbin1;
Fig=figure(2*nt+10);
pp1=histogram(x_sample(:,2,nt),nbin2,'Normalization','pdf');
hold on
pp2=plot(x2,ypdf_ini(:,2),'r-.','LineWidth',2.0);
hold on
pp3=plot(x2,ypdf_cal(:,2),'k-','LineWidth',2.0);
% hold on
% pp3=plot(x2,ypdf_mle(:,2),'b:','LineWidth',2.0);
xlim([min(CI(2,1),CI_(2,1))-1000, max(CI(2,2),CI_(2,2))-1000])
ylim([0 7]*10^-4)
set(gca,'fontsize', 13,'fontname','times');
ytickangle(90)
xlabel('Density (kg/m^{3}) - \theta_{2}','fontsize', 15)
ylabel('PDF','fontsize', 15)
% legend([pp1 pp2 pp3],'Virtual Input Data','Prior PDF','Calibrated PDF','Location','best','fontsize', 12)
box on
set(Fig,'pos',[100 541 660 303]);
    end
end