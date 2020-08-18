clc
close all
clear

%% Mathematical Model (Corr)
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

del_size=[1:1:20];
count=[1:size(del_size,2)];
x_axis_log=10.^(-[del_size(count)]);
plot_axis=[1:size(del_size,2)];
for nt=1:2
    file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
        sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200229.mat'],dpt_n(2),dpt_n(1),1)];
    [ff_analy, gg_analy]=MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(x0,file_exp_name,file_s_name,type_num(nt,:),perf_func); % Likelihood Function Metric(KDE)
    for ii=count
        [ff_cvm(ii,:), gg_cvm(ii,:)]=MLE_kerkrig_Ndim_numer_cvm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt,:),perf_func);
        [ff_fdm(ii,:), gg_fdm(ii,:)]=MLE_kerkrig_Ndim_numer_fdm_corr_parfor(x0,10^-del_size(ii),file_exp_name,file_s_name,type_num(nt,:),perf_func);
        disp(sprintf('gradient calculation end %d',ii))
    end
    err_cvm(:,:,nt)=abs((gg_cvm-gg_analy)./gg_analy)*100;
    err_fdm(:,:,nt)=abs((gg_fdm-gg_analy)./gg_analy)*100;
    disp(sprintf('%d번째 케이스',nt));
end

save('SA_vs_CVM_corr_200427.mat','err_cvm','err_fdm','gg_cvm','gg_analy','ff_cvm','ff_analy','gg_fdm','ff_fdm')

for kk=1%:4
    Fig=figure(kk);
    set(Fig,'pos',[1 541 820 493]);
    p1=loglog(x_axis_log(count),err_fdm(plot_axis(count),1,kk),'r-o','LineWidth',1.5);
    hold on
    p2=loglog(x_axis_log(count),err_fdm(plot_axis(count),2,kk),'r:o','LineWidth',1.5);
    hold on
    p3=loglog(x_axis_log(count),err_fdm(plot_axis(count),3,kk),'r-.o','LineWidth',1.5);
    hold on
    p4=loglog(x_axis_log(count),err_fdm(plot_axis(count),4,kk),'r--o','LineWidth',1.5);
    hold on
    p5=loglog(x_axis_log(count),err_cvm(plot_axis(count),1,kk),'k-d','LineWidth',1.5);
    hold on
    p6=loglog(x_axis_log(count),err_cvm(plot_axis(count),2,kk),'k:d','LineWidth',1.5);
    hold on
    p7=loglog(x_axis_log(count),err_cvm(plot_axis(count),3,kk),'k-.d','LineWidth',1.5);
    hold on
    p8=loglog(x_axis_log(count),err_cvm(plot_axis(count),4,kk),'k--d','LineWidth',1.5);
    box on
    grid on
    xlim([10^-20 10^-1])
    ylim([10^-14 10^4])
    axx = gca;
    axx.FontSize = 12;
    xlabel('Perturbation','FontSize',15)
    ylabel('Error (%)','FontSize',15)
    hleg1=legend([p1,p2,p3,p4,p5,p6,p7,p8], ...
        {'FDM - \it$\frac{df}{d\mu_{1}}$','FDM - \it$\frac{df}{d\sigma_{1}}$','FDM - \it$\frac{df}{d\mu_{2}}$','FDM - \it$\frac{df}{d\sigma_{2}}$', ...
        'CVM - \it$\frac{df}{d\mu_{1}}$','CVM - \it$\frac{df}{d\sigma_{1}}$','CVM - \it$\frac{df}{d\mu_{2}}$','CVM - \it$\frac{df}{d\sigma_{2}}$'},'location','northeastoutside','FontSize',18,'Interpreter','latex');
    set(gca,'fontname','times')  % Set it to times
end

err_table_fdm(:,:)=err_fdm(8,:,:);
err_table_cvm(:,:)=err_cvm(15,:,:);
err_table_fdm=err_table_fdm';
err_table_cvm=err_table_cvm';
