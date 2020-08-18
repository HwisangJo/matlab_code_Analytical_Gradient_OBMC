% 20200221 : 모르는 입력변수의 신뢰구간에서 (1) Max MSE, (2) Adaptive 방법 중 어느 것이 수렴 속도가 더 빠른지 검증
% model_ex_ver02, 아웃풋 1개, initial samples 7개, 더 넓은 범위
% (2)
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
%%
foldername5=[opt_folder,'\','opt_hist_analy'];
mkdir(foldername5)
for nt=2%:2
load([pwd,'\opt_data_ex2_krig_200309\opt_input_mle.mat'])
    folder_name='\exp_data_ex2_bracket_200316';
% poolobj=parpool(8);
    for ini=1:2%:5
for i=1:6

% foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d.mat', nt)];
% foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_norm.mat', nt)];
foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_norm_trr_analy.mat', nt,ini,i)];
% foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_norm_sqp.mat', nt)];

% Optimization-Based Statistical Model Calibration (Parameter Estimation)
file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
    sprintf(['exp_data_ex2_krig_',dist_name_type{nt},'_n%d_%d_trial%d_200309.mat'],dpt_n(2),dpt_n(1),1)];
file_s_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',i)];
save_sign=0;

% Performance Function
perf_func1=@(x) model_krig_modal_bracket(x,2);
perf_func2=@(x) model_krig_modal_bracket2(x,1);
perf_func=@(x,y_num) model_krig_mult_resp(x,y_num,perf_func1,perf_func2);
% x_true=[210e9, 210e9*0.05, 8000, 8000*0.055];
% lbx=[2.5e7, 2.8e7*0.045, 6500, 8000*0.045];
% ubx=[3.1e7, 2.8e7*0.055, 9500, 8000*0.055];

% lbx=[210e9, 210e9*0.05, 8000, 8000*0.05]*0.5;
% ubx=[210e9, 210e9*0.05, 8000, 8000*0.05]*1.5;

% lbx=[210e9, 210e9*0.05, 8000, 8000*0.1]*0.5;
% ubx=[210e9, 210e9*0.05, 8000, 8000*0.1]*1.5;

lbx=[210e9, 210e9*0.10, 8000, 8000*0.10]*0.5;
ubx=[210e9, 210e9*0.10, 8000, 8000*0.10]*1.5;

% x0=(lbx+ubx)/2*0.9;
ratio_x0=[1.10 1.05 0.90 0.95;0.92 0.98 1.10 1.10; 0.95 0.95 1.10 1.05;1.03 1.07 0.98 0.92;1.04 0.95 0.92 0.90];
x0=parm_true2(nt,:).*ratio_x0(ini,:);

% x0=[2.7e7, 210e9*0.08, 8500, 8000*0.08];
% x0=[3.0e7, 2.8e7*0.052, 8000, 8000*0.052];
% x0=[210e9, 1.3e6, 8.0e3, 390];
% x0_tot=[210e9, 1.3e6, 8.0e3, 390;
%     210e9, 1.3e6, 8.1e3, 350];
% x0_tot=[210e9, 1.3e6, 8.2e3, 810;
%     210e9, 1.3e6, 8.1e3, 650];
% abs(x0-x_true)./x_true*100
% x0=x0_tot(nt,:);

lbx_n=zeros(1,size(lbx,2));
ubx_n=ones(1,size(ubx,2));
% x0_n=0.5*ones(1,size(ubx,2));
x0_n=(x0-lbx)./(ubx-lbx);

% obj = @(x) MLE_kerkrig_Ndim_analy_mathe(x,file_exp_name,file_s_name,type_num(nt),perf_func); % Likelihood Function Metric(KDE)
% obj = @(x) MLE_kerkrig_Ndim_analy_mathe_parfor(x,file_exp_name,file_s_name,type_num(nt),perf_func); % Likelihood Function Metric(KDE)
% obj = @(x) MLE_kerkrig_Ndim_analy_krig(x,file_exp_name,file_s_name,type_num(nt),perf_func,save_sign); % Likelihood Function Metric(KDE)
% obj = @(x) MLE_kerkrig_Ndim_analy_krig_parfor(x,file_exp_name,file_s_name,type_num_case(nt,:),perf_func,save_sign); % Likelihood Function Metric(KDE)
obj = @(x) MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(x,file_exp_name,file_s_name,type_num_case(nt,:),perf_func); % Likelihood Function Metric(KDE)
obj_n = @(x) perf_norm(x,obj,lbx,ubx);

% Optimization
tic
opts = optimoptions('fmincon','Display','iter','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'StepTolerance',1e-15,'OptimalityTolerance',1e-15,'ConstraintTolerance',1e-15,'MaxFunctionEvaluations',1000);
% opts = optimoptions('fmincon','Display','iter','Algorithm','sqp','SpecifyObjectiveGradient',true,'StepTolerance',1e-15,'OptimalityTolerance',1e-15,'ConstraintTolerance',1e-15,'MaxFunctionEvaluations',1000);
% opts = optimoptions('fmincon','Display','iter','Algorithm','sqp','SpecifyObjectiveGradient',true,'CheckGradients',true,'FiniteDifferenceType','central','StepTolerance',1e-15);
% [xvar,fval_sqp,exitflag,output] = fmincon(obj,x0,[],[],[],[],lbx,ubx,[],opts);
% [xvar_n,fval_sqp,exitflag,output] = fmincon(obj_n,x0_n,[],[],[],[],lbx_n,ubx_n,[],opts);
[xvar_n,fval_sqp,history] = runfmincon(obj_n,x0_n,lbx_n,ubx_n);
% [xvar_n,fval_sqp,history] = myproblem(obj_n,x0_n,lbx_n,ubx_n,opts);
xvar=xvar_n.*(ubx-lbx)+lbx;
XXX=[xvar(1) xvar(2) xvar(3) xvar(4) fval_sqp]
time_opt_sqp(i,:)=toc

save(foldername3,'xvar','fval_sqp','time_opt_sqp','lbx','ubx','history','x0') % sim?

% delete(gcp('nocreate'))

load([pwd,'\opt_data_ex2_krig_200309\opt_input_mle.mat'])
% scatter(history.x())
x_un=history.x.*(ubx-lbx)+lbx;
% figure

Fig=figure(nt);
set(Fig,'pos',[100 541 2260 743]);
ax=subplot(1,2,1);
scatter(x_un(1,1),x_un(1,2),'ro')
hold on
scatter(x_un(end,1),x_un(end,2),'rx')
hold on
plot(x_un(:,1),x_un(:,2),'k.-')
box on
hold on
scatter(parm_true2(nt,1),parm_true2(nt,2),'bd')
xlabel('Mean of \theta_{1}')
ylabel('STD of \theta_{1}')
% title('History of SQP - \theta_{1}')
title('History of TRR - \theta_{1}')
xlim([lbx(1) ubx(1)])
ylim([lbx(2) ubx(2)])

ax=subplot(1,2,2);
scatter(x_un(1,3),x_un(1,4),'ro')
hold on
scatter(x_un(end,3),x_un(end,4),'rx')
hold on
plot(x_un(:,3),x_un(:,4),'k.-')
box on
hold on
scatter(parm_true2(nt,3),parm_true2(nt,4),'bd')
xlabel('Mean of \theta_{2}')
ylabel('STD of \theta_{2}')
% title('History of SQP - \theta_{2}')
title('History of TRR - \theta_{2}')
xlim([lbx(3) ubx(3)])
ylim([lbx(4) ubx(4)])
saveas(Fig,[foldername5,'\',sprintf('optimum_history_case%d_initial%d_seed%d',nt,ini,i),'.jpg'])
close all
end
    end
end


















