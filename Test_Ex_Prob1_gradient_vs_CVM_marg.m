clc
close all
clear

%% Mathematical Model (Marginal)
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

% Initial Points
x0=[16.0 1.5];

% Performance Function
perf_func=@(x) model_ex_1d_nonlinear(x);

del_size=[1:1:20];
count=[1:size(del_size,2)];
x_axis_log=10.^(-[del_size(count)]);
plot_axis=[1:size(del_size,2)];
for nt=1:4
    file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
        sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(2),dpt_n(1),1)];
    [ff_analy, gg_analy]=MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(x0,file_exp_name,file_s_name,type_num(nt),perf_func); % Likelihood Function Metric(KDE)
    for ii=count
        [ff_cvm(ii,:), gg_cvm(ii,:)]=MLE_kerkrig_Ndim_numer_cvm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt),perf_func);
        [ff_fdm(ii,:), gg_fdm(ii,:)]=MLE_kerkrig_Ndim_numer_fdm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt),perf_func);
        disp(sprintf('gradient calculation end %d',ii))
    end
    err_cvm(:,:,nt)=abs((gg_cvm-gg_analy)./gg_analy)*100;
    err_fdm(:,:,nt)=abs((gg_fdm-gg_analy)./gg_analy)*100;
    disp(sprintf('%d th case',nt));
end

save('SA_vs_CVM_marg_200427.mat','err_cvm','err_fdm','gg_cvm','gg_analy','ff_cvm','ff_analy','gg_fdm','ff_fdm')

for kk=1%:4
    Fig=figure(kk);
    set(Fig,'pos',[1 541 810 493]);
    p1=loglog(x_axis_log(count),err_fdm(plot_axis(count),1,kk),'r-o','LineWidth',1.5);
    hold on
    p2=loglog(x_axis_log(count),err_fdm(plot_axis(count),2,kk),'r:o','LineWidth',1.5);
    hold on
    p5=loglog(x_axis_log(count),err_cvm(plot_axis(count),1,kk),'k-d','LineWidth',1.5);
    hold on
    p6=loglog(x_axis_log(count),err_cvm(plot_axis(count),2,kk),'k:d','LineWidth',1.5);
    box on
    grid on
    xlim([10^-20 10^-1])
    ylim([10^-14 10^4])
    axx = gca;
    axx.FontSize = 12;
    xlabel('Perturbation','FontSize',15)
    ylabel('Error (%)','FontSize',15)
    legend([p1,p2,p5,p6], ...
        {'FDM - \it$\frac{df}{d\mu}$','FDM - \it$\frac{df}{d\sigma}$', ...
        'CVM - \it$\frac{df}{d\mu}$','CVM - \it$\frac{df}{d\sigma}$'},'location','northeastoutside','FontSize',18,'Interpreter','latex')
    set(gca,'fontname','times')  % Set it to times
end

err_table_fdm(:,:)=err_fdm(8,:,:);
err_table_cvm(:,:)=err_cvm(15,:,:);
err_table_fdm=err_table_fdm';
err_table_cvm=err_table_cvm';











