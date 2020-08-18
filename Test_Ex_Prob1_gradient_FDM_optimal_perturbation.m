clc
close all
clear

%% Mathematical Model (Univariate-linear)
file_s_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',1)];
folder_name_type={'\exp_data_case1_normal_krig_200222';
    '\exp_data_case2_lognormal_krig_200222';
    '\exp_data_case3_gumbel_krig_200222';
    '\exp_data_case4_extreme_krig_200222';
    };
dist_type={'Normal','Lognormal','Extreme Value','Extreme Value'};
dist_name_type={'case1_normal';'case2_lognormal';'case3_gumbel';'case4_extreme';'case5_uniform';};
type_num=[1 2 3 6 5];
dpt_n=[50 1;]; % dimension and number of experiment results

% folder_name='\exp_data_ex1_nonlinear_200316';
folder_name='\exp_data_ex1_linear_200316';

% Initial Points
x0=[16.0 1.5];

% Performance Function
% perf_func=@(x) model_ex_1d_nonlinear(x);
perf_func=@(x) model_ex_1d_linear(x);

% poolobj=parpool(8);
del_size=[1:1:20];
count=[1:size(del_size,2)];
x_axis_log=10.^(-[del_size(count)]);
plot_axis=[1:size(del_size,2)];
for nt=1:4
file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
    sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(2),dpt_n(1),1)];
    [ff_analy, gg_analy]=MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(x0,file_exp_name,file_s_name,type_num(nt),perf_func); % Likelihood Function Metric(KDE)
    for ii=count
%         [ff_cvm(ii,:), gg_cvm(ii,:)]=MLE_kerkrig_Ndim_numer_cvm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt),perf_func);
        [ff_fdm(ii,:), gg_fdm(ii,:)]=MLE_kerkrig_Ndim_numer_fdm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt),perf_func);
        disp(sprintf('gradient calculation end %d',ii))
    end
%     err_cvm(:,:,nt)=abs((gg_cvm-gg_analy)./gg_analy)*100;
    err_fdm(:,:,nt)=abs((gg_fdm-gg_analy)./gg_analy)*100;
    disp(sprintf('%d번째 케이스',nt));
end
% delete(gcp('nocreate'))

Err_fdm=mean(err_fdm,2);
for kk=1:4
Fig=figure;
set(Fig,'pos',[1 541 810 493]);
p1=loglog(x_axis_log(count),err_fdm(plot_axis(count),1,kk),'r-o','LineWidth',1.5);
hold on
p2=loglog(x_axis_log(count),err_fdm(plot_axis(count),2,kk),'r:o','LineWidth',1.5);
hold on
p3=loglog(x_axis_log(count),Err_fdm(plot_axis(count),1,kk),'k-o','LineWidth',1.5);
% hold on
% p5=loglog(x_axis_log(count),err_cvm(plot_axis(count),1,kk),'k-d','LineWidth',1.5);
% hold on
% p6=loglog(x_axis_log(count),err_cvm(plot_axis(count),2,kk),'k:d','LineWidth',1.5);
box on
grid on
xlim([10^-20 10^-1])
% ylim([10^-14 10^4])
axx = gca;
axx.FontSize = 12;
xlabel('Perturbation','FontSize',15)
ylabel('Error (%)','FontSize',15)
legend([p1,p2], ... 
    {'FDM - df/d\mu','FDM - df/d\sigma' ...
    },'location','northeastoutside','FontSize',12)
set(gca,'fontname','times')  % Set it to times
end

%% Mathematical Model (Univariate-nonlinear)
file_s_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',1)];
folder_name_type={'\exp_data_case1_normal_krig_200222';
    '\exp_data_case2_lognormal_krig_200222';
    '\exp_data_case3_gumbel_krig_200222';
    '\exp_data_case4_extreme_krig_200222';
    };
dist_type={'Normal','Lognormal','Extreme Value','Extreme Value'};
dist_name_type={'case1_normal';'case2_lognormal';'case3_gumbel';'case4_extreme';'case5_uniform';};
type_num=[1 2 3 6 5];
dpt_n=[50 1;]; % dimension and number of experiment results

folder_name='\exp_data_ex1_nonlinear_200316';
% folder_name='\exp_data_ex1_linear_200316';

% Initial Points
x0=[16.0 1.5];

% Performance Function
perf_func=@(x) model_ex_1d_nonlinear(x);
% perf_func=@(x) model_ex_1d_linear(x);

% poolobj=parpool(8);
del_size=[1:1:20];
count=[1:size(del_size,2)];
x_axis_log=10.^(-[del_size(count)]);
plot_axis=[1:size(del_size,2)];
for nt=1:4
file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
    sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(2),dpt_n(1),1)];
    [ff_analy, gg_analy]=MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(x0,file_exp_name,file_s_name,type_num(nt),perf_func); % Likelihood Function Metric(KDE)
    for ii=count
%         [ff_cvm(ii,:), gg_cvm(ii,:)]=MLE_kerkrig_Ndim_numer_cvm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt),perf_func);
        [ff_fdm(ii,:), gg_fdm(ii,:)]=MLE_kerkrig_Ndim_numer_fdm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt),perf_func);
        disp(sprintf('gradient calculation end %d',ii))
    end
%     err_cvm(:,:,nt)=abs((gg_cvm-gg_analy)./gg_analy)*100;
    err_fdm(:,:,nt)=abs((gg_fdm-gg_analy)./gg_analy)*100;
    disp(sprintf('%d번째 케이스',nt));
end
% delete(gcp('nocreate'))

Err_fdm=mean(err_fdm,2);
for kk=1:4
Fig=figure;
set(Fig,'pos',[1 541 810 493]);
p1=loglog(x_axis_log(count),err_fdm(plot_axis(count),1,kk),'r-o','LineWidth',1.5);
hold on
p2=loglog(x_axis_log(count),err_fdm(plot_axis(count),2,kk),'r:o','LineWidth',1.5);
hold on
p3=loglog(x_axis_log(count),Err_fdm(plot_axis(count),1,kk),'k-o','LineWidth',1.5);
% hold on
% p5=loglog(x_axis_log(count),err_cvm(plot_axis(count),1,kk),'k-d','LineWidth',1.5);
% hold on
% p6=loglog(x_axis_log(count),err_cvm(plot_axis(count),2,kk),'k:d','LineWidth',1.5);
box on
grid on
xlim([10^-20 10^-1])
% ylim([10^-14 10^4])
axx = gca;
axx.FontSize = 12;
xlabel('Perturbation','FontSize',15)
ylabel('Error (%)','FontSize',15)
legend([p1,p2], ... 
    {'FDM - df/d\mu','FDM - df/d\sigma' ...
    },'location','northeastoutside','FontSize',12)
set(gca,'fontname','times')  % Set it to times
end

%% Mathematical Model (Bivariate)
dpt_n=[50 1]; % dimension and number of experiment results
file_s_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',1)];
folder_name_type={'\exp_data_case1_normal_lognormal_200229';
    '\exp_data_case2_gumbel_extreme_200229';
    };
dist_type={'Normal','Lognormal';'Extreme Value','Extreme Value'};
dist_name_type={'case1_normal_lognormal';'case2_gumbel_extreme';};
type_num=[1 2; 3 6];
folder_name='\exp_data_ex1_exp_multi_200316';

% Initial Points
x0=[11.0 1.0 11.0 1.0];

% Performance Function
perf_func=@(x,y_num) model_ex_mult_resp(x,y_num);

% poolobj=parpool(8);
del_size=[1:1:20];
count=[1:size(del_size,2)];
x_axis_log=10.^(-[del_size(count)]);
plot_axis=[1:size(del_size,2)];
for nt=1:2
    file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
        sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200229.mat'],dpt_n(2),dpt_n(1),1)];
    [ff_analy, gg_analy]=MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(x0,file_exp_name,file_s_name,type_num(nt,:),perf_func); % Likelihood Function Metric(KDE)
    for ii=count
%         [ff_cvm(ii,:), gg_cvm(ii,:)]=MLE_kerkrig_Ndim_numer_cvm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt,:),perf_func);
        [ff_fdm(ii,:), gg_fdm(ii,:)]=MLE_kerkrig_Ndim_numer_fdm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt,:),perf_func);
        disp(sprintf('gradient calculation end %d',ii))
    end
%     err_cvm(:,:,nt)=abs((gg_cvm-gg_analy)./gg_analy)*100;
    err_fdm(:,:,nt)=abs((gg_fdm-gg_analy)./gg_analy)*100;
    disp(sprintf('%d번째 케이스',nt));
end
% delete(gcp('nocreate'))

Err_fdm=mean(err_fdm,2);
for kk=1:2
Fig=figure;
set(Fig,'pos',[1 541 810 493]);
p1=loglog(x_axis_log(count),err_fdm(plot_axis(count),1,kk),'r-o','LineWidth',1.5);
hold on
p2=loglog(x_axis_log(count),err_fdm(plot_axis(count),2,kk),'r:o','LineWidth',1.5);
hold on
p3=loglog(x_axis_log(count),err_fdm(plot_axis(count),3,kk),'r-.o','LineWidth',1.5);
hold on
p4=loglog(x_axis_log(count),err_fdm(plot_axis(count),4,kk),'r--o','LineWidth',1.5);
hold on
p5=loglog(x_axis_log(count),Err_fdm(plot_axis(count),1,kk),'k-o','LineWidth',1.5);
% hold on
% p5=loglog(x_axis_log(count),err_cvm(plot_axis(count),1,kk),'k-d','LineWidth',1.5);
% hold on
% p6=loglog(x_axis_log(count),err_cvm(plot_axis(count),2,kk),'k:d','LineWidth',1.5);
% hold on
% p7=loglog(x_axis_log(count),err_cvm(plot_axis(count),3,kk),'k-.d','LineWidth',1.5);
% hold on
% p8=loglog(x_axis_log(count),err_cvm(plot_axis(count),4,kk),'k--d','LineWidth',1.5);
box on
grid on
xlim([10^-20 10^-1])
% ylim([10^-14 10^4])
axx = gca;
axx.FontSize = 12;
xlabel('Perturbation','FontSize',15)
ylabel('Error (%)','FontSize',15)
legend([p1,p2,p3,p4], ... 
    {'FDM - df/d\mu_{1}','FDM - df/d\sigma_{1}','FDM - df/d\mu_{2}','FDM - df/d\sigma_{2}' ...
    },'location','northeastoutside','FontSize',12)
set(gca,'fontname','times')  % Set it to times
end

%% Bracket Model (Univariate)
addpath(genpath([pwd,'\imm1460\dace']))
dpt_n=[50 1]; % dimension and number of experiment results
file_s_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',1)];
folder_name_type={'\exp_data_case1_normal_krig_200222';
    '\exp_data_case2_lognormal_krig_200222';
    '\exp_data_case3_gumbel_krig_200222';
    '\exp_data_case4_extreme_krig_200222';
    };
dist_type={'Normal','Lognormal','Extreme Value','Extreme Value'};
dist_name_type={'case1_normal';'case2_lognormal';'case3_gumbel';'case4_extreme';'case5_uniform';};
type_num=[1 2 3 6 5];
dpt_n=[50 1;]; % dimension and number of experiment results
folder_name='\exp_data_ex2_bracket_1d_200316';

% Initial Points
x0=[210e9, 210e9*0.1].*[1.1 1.05];

% Performance Function
perf_func=@(x) model_krig_modal_bracket2(x,1);

% poolobj=parpool(8);
del_size=[1:1:20];
count=[1:size(del_size,2)];
x_axis_log=10.^(-[del_size(count)]);
plot_axis=[1:size(del_size,2)];
for nt=2:4
for trial=1
    file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
        sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(2),dpt_n(1),1)];
    [ff_analy, gg_analy]=MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(x0,file_exp_name,file_s_name,type_num,perf_func);
    for ii=count
%     time_grad_cvm(nt,trial)=0;gg_cvm=0;
%     [~,gg_cvm,time_grad_cvm(nt,trial)]=MLE_kerkrig_Ndim_numer_cvm_corr_parfor(x0,10^-del_size(15),file_exp_name,file_s_name,type_num(nt,:),perf_func);
    [ff_fdm(ii,:), gg_fdm(ii,:)]=MLE_kerkrig_Ndim_numer_fdm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num,perf_func);
    disp(['%%%%%%%%%%%%%%%%',sprintf(' %d번째 perturbation ',ii),'%%%%%%%%%%%%%%%%']);
    end
    err_fdm(:,:,nt)=abs((gg_fdm-gg_analy)./gg_analy)*100;
    
% disp(['%%%%%%%%%%%%%%%%',sprintf(' %d번째 트라이얼 ',trial),'%%%%%%%%%%%%%%%%']);
end
disp(['%%%%%%%%%%%%%%%%',sprintf(' %d번째 분포 ',nt),'%%%%%%%%%%%%%%%%']);
end
% delete(gcp('nocreate'))

Err_fdm=mean(err_fdm,2);
for kk=1:4
figure(kk)
p1=loglog(x_axis_log(count),err_fdm(plot_axis(count),1,kk),'r-.d','LineWidth',1.5);
hold on
p2=loglog(x_axis_log(count),err_fdm(plot_axis(count),2,kk),'r:o','LineWidth',1.5);
hold on
p3=loglog(x_axis_log(count),Err_fdm(plot_axis(count),1,kk),'k-x','LineWidth',1.5);
box on
grid on
% xlim([10^-20 10^-1])
% ylim([10^-14 10^4])
axx = gca;
axx.FontSize = 12;
xlabel('Perturbation','FontSize',15)
ylabel('Error (%)','FontSize',15)
legend([p1,p2],{'Error of df/d\mu','Error of df/d\sigma'},'location','best','FontSize',12)
end

%% Bracket Model (Bivariate)
addpath(genpath([pwd,'\imm1460\dace']))

dpt_n=[50 1]; % dimension and number of experiment results
file_s_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',1)];
folder_name_type={'\exp_data_case1_normal_lognormal_krig_200309';
    '\exp_data_case2_gumbel_extreme_krig_200309';
    };
dist_type={'Normal','Lognormal';'Extreme Value','Extreme Value'};
dist_name_type={'case1_normal_lognormal';'case2_gumbel_extreme';};
type_num=[1 2; 3 6];
folder_name='\exp_data_ex2_bracket_200316';

% Initial Points (Predetermined Point)
x0=[210e9, 210e9*0.1, 8000, 8000*0.1].*[1.1 1.05 0.9 0.95];

% Performance Function
perf_func1=@(x) model_krig_modal_bracket(x,2);
perf_func2=@(x) model_krig_modal_bracket2(x,1);
perf_func=@(x,y_num) model_krig_mult_resp(x,y_num,perf_func1,perf_func2);

% poolobj=parpool(8);
del_size=[1:1:20];
count=[1:size(del_size,2)];
x_axis_log=10.^(-[del_size(count)]);
plot_axis=[1:size(del_size,2)];
for nt=2%:2
for trial=1%:10
    file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
        sprintf(['exp_data_ex2_krig_',dist_name_type{nt},'_n%d_%d_trial%d_200309.mat'],dpt_n(2),dpt_n(1),1)];
    [ff_analy, gg_analy]=MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(x0,file_exp_name,file_s_name,type_num(nt,:),perf_func);
    for ii=count
    time_grad_cvm(nt,trial)=0;gg_cvm=0;
%     [~,gg_cvm,time_grad_cvm(nt,trial)]=MLE_kerkrig_Ndim_numer_cvm_corr_parfor(x0,10^-del_size(15),file_exp_name,file_s_name,type_num(nt,:),perf_func);
    [ff_fdm(ii,:), gg_fdm(ii,:)]=MLE_kerkrig_Ndim_numer_fdm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt,:),perf_func);
    disp(['%%%%%%%%%%%%%%%%',sprintf(' %d번째 perturbation ',ii),'%%%%%%%%%%%%%%%%']);
    end
%     err_cvm{trial,nt}=abs((gg_cvm-gg_analy)./gg_analy)*100;
%     err_fdm{trial,nt}=abs((gg_fdm-gg_analy)./gg_analy)*100;
    err_fdm(:,:,nt)=abs((gg_fdm-gg_analy)./gg_analy)*100;
%     clear gg gg_
% disp(['%%%%%%%%%%%%%%%%',sprintf(' %d번째 트라이얼 ',trial),'%%%%%%%%%%%%%%%%']);
end
disp(['%%%%%%%%%%%%%%%%',sprintf(' %d번째 분포 ',nt),'%%%%%%%%%%%%%%%%']);
end
% delete(gcp('nocreate'))

Err_fdm=mean(err_fdm,2);
for kk=2%:2
Fig=figure;
set(Fig,'pos',[1 541 810 493]);
p1=loglog(x_axis_log(count),err_fdm(plot_axis(count),1,kk),'r-o','LineWidth',1.5);
hold on
p2=loglog(x_axis_log(count),err_fdm(plot_axis(count),2,kk),'r:o','LineWidth',1.5);
hold on
p3=loglog(x_axis_log(count),err_fdm(plot_axis(count),3,kk),'r-.o','LineWidth',1.5);
hold on
p4=loglog(x_axis_log(count),err_fdm(plot_axis(count),4,kk),'r--o','LineWidth',1.5);
hold on
p5=loglog(x_axis_log(count),Err_fdm(plot_axis(count),1,kk),'k-x','LineWidth',1.5);
% hold on
% p5=loglog(x_axis_log(count),err_cvm(plot_axis(count),1,kk),'k-d','LineWidth',1.5);
% hold on
% p6=loglog(x_axis_log(count),err_cvm(plot_axis(count),2,kk),'k:d','LineWidth',1.5);
% hold on
% p7=loglog(x_axis_log(count),err_cvm(plot_axis(count),3,kk),'k-.d','LineWidth',1.5);
% hold on
% p8=loglog(x_axis_log(count),err_cvm(plot_axis(count),4,kk),'k--d','LineWidth',1.5);
box on
grid on
xlim([10^-20 10^-1])
% ylim([10^-14 10^4])
axx = gca;
axx.FontSize = 12;
xlabel('Perturbation','FontSize',15)
ylabel('Error (%)','FontSize',15)
% legend([p1,p2,p3,p4,p5,p6,p7,p8], ... 
%     {'Error of df/d\mu_{1} - FDM','Error of df/d\sigma_{1} - FDM','Error of df/d\mu_{2} - FDM','Error of df/d\sigma_{2} - FDM', ...
%     'Error of df/d\mu_{1} - CVM','Error of df/d\sigma_{1} - CVM','Error of df/d\mu_{2} - CVM','Error of df/d\sigma_{2} - CVM'},'location','northeastoutside','FontSize',12)
legend([p1,p2,p3,p4,p5], ... 
    {'FDM - df/d\mu_{1}','FDM - df/d\sigma_{1}','FDM - df/d\mu_{2}','FDM - df/d\sigma_{2}','FDM - Avg' ...
    },'location','northeastoutside','FontSize',12)
set(gca,'fontname','times')  % Set it to times
end


