clc
close all
clear

addpath(genpath([pwd,'\imm1460\dace']))
%% Mathematical Model
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
x0=[210e9, 210e9*0.05];

% Performance Function
perf_func=@(x) model_krig_modal_bracket2(x,1);

del_size=[1:1:20];
count=[1:size(del_size,2)];
x_axis_log=10.^(-[del_size(count)]);
plot_axis=[1:size(del_size,2)];
for nt=1%:4
    for trial=1:50
        file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
            sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(2),dpt_n(1),1)];
        
        [~,gg_cvm,time_grad_cvm(nt,trial)]=MLE_kerkrig_Ndim_numer_cvm_corr_parfor(x0,10^-del_size(15),file_exp_name,file_s_name,type_num(nt),perf_func);
        [~,gg_fdm,time_grad_fdm(nt,trial)]=MLE_kerkrig_Ndim_numer_fdm_corr_parfor(x0,10^-del_size(7),file_exp_name,file_s_name,type_num(nt),perf_func);
        [~,gg_analy,time_grad_analy(nt,trial)]=MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(x0,file_exp_name,file_s_name,type_num(nt),perf_func);
        
        err_cvm{trial,nt}=abs((gg_cvm-gg_analy)./gg_analy)*100;
        err_fdm{trial,nt}=abs((gg_fdm-gg_analy)./gg_analy)*100;
        
        disp(['%%%%%%%%%%%%%%%%',sprintf(' %d th trial ',trial),'%%%%%%%%%%%%%%%%']);
    end
    disp(['%%%%%%%%%%%%%%%%',sprintf(' %d th distribution ',nt),'%%%%%%%%%%%%%%%%']);
end

time_total=[mean(time_grad_fdm,2) mean(time_grad_cvm,2) mean(time_grad_analy,2)];

% save('time_gradient_vs_CVM_FDM_bracket_1d_200427_10_3.mat')
% save('time_gradient_vs_CVM_FDM_bracket_1d_200427_10_4.mat')
% save('time_gradient_vs_CVM_FDM_bracket_1d_200427_10_5.mat')
save('time_gradient_vs_CVM_FDM_bracket_1d_200427_10_6.mat')

