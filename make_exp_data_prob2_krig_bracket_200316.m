% 2020.02.22 Make Experiment Data
%%
clc
close all
clear

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
%% 가상 시험 데이터 생성, 2020.02.22, 범위 넓은 조건
dpt_n=[50 1;]; % dimension and number of experiment results
% xc_mean_(1)=210e9;xc_std_(1)=0.05*xc_mean_(1);
% xc_mean_(2)=8000;xc_std_(2)=0.05*xc_mean_(2);
xc_mean_(1)=210e9;xc_std_(1)=0.10*xc_mean_(1);
xc_mean_(2)=8000;xc_std_(2)=0.10*xc_mean_(2);

perf_func1=@(x) model_krig_modal_bracket(x,2);
perf_func2=@(x) model_krig_modal_bracket2(x,1);
for nt=1%:2
    folder_name='\exp_data_ex2_bracket_200316';
    mkdir([pwd,folder_name,folder_name_type{nt}])
    nexp=1;
for trial=1:nexp
    for pt=1
        % Observation Sites
        ne=dpt_n(pt,2);tr=dpt_n(pt,1); % dimension and number of experiment results     
        for ui=1:2
        ab(ui,:)=convert_param(xc_mean_(ui),xc_std_(ui),type_num_case(nt,ui));
        xc_mean(ui)=ab(ui,1);xc_std(ui)=ab(ui,2);
        
        if type_num_case(nt,ui)==3
            xc(:,ui,:)=xc_mean(ui)-(random(dist_type{nt,ui},xc_mean(ui),xc_std(ui), [ne,tr])-xc_mean(ui)); % Unknown Variables
        else
            xc(:,ui,:)=random(dist_type{nt,ui},xc_mean(ui),xc_std(ui), [ne,tr]); % Unknown Variables
        end
        end
        
        % Model
        y_exp=[];x_exp=[];x_v=[];
        for ii=1:2
        for i=1:tr
            v_exp_data(:,ii,i)=model_krig_mult_resp([xc(:,:,i)],ii,perf_func1,perf_func2);
            x_v(i,ii)=xc(:,ii,i);
            y_exp=[y_exp;v_exp_data(:,ii,i)];
        end
        end
        
        % Sample mena, variance
        mean_samp(1,:)=mean(x_v);
        sd_samp(1,:)=std(x_v);
        var_samp(1,:)=var(x_v);
        
        % 저장
        file_name=[pwd,folder_name,folder_name_type{nt},'\', ...
            sprintf(['exp_data_ex2_krig_',dist_name_type{nt},'_n%d_%d_trial%d_200309.mat'],dpt_n(pt,2),dpt_n(pt,1),trial)];
%         save(file_name,'v_exp_data','xc','ne','tr','y_exp','mean_samp','sd_samp','var_samp','x_v') %
%         clear v_exp_data xc ne tr y_exp x_v
    end
end
end

% lbx=[7, 10^-3];
% ubx=[13, 2.5];
% x0=(lbx+ubx)/2;
% file_s_name=[pwd,'\s_save_191017','\',sprintf('s_save_trial%d_191017.mat',1)];
% % perf_func=@(x) model_ex_ver02(x);
% % perf_func=@(x) branin(x);
% % perf_func=@(x) model_krig(x);
% perf_func=@(x) model_krig_initial(x);
% % ff=MLE_kerkrig_Ndim_analy_mathe(x0,file_name,file_s_name,1,perf_func);
% % ff=MLE_kerkrig_Ndim_analy_krig(x0,file_name,file_s_name,1,perf_func,0);

xc_mean_(1)=2.8e7;xc_std_(1)=0.05*xc_mean_(1);
xc_mean_(2)=7800;xc_std_(2)=0.05*xc_mean_(2);

figure(2)
scatter(repmat(xd,tr,1),y_exp(:,1),'k.')
axis([0.06 0.09 2.5 7])
box on

figure(3)
scatter(x_v(:,1),x_v(:,2),'k.')
CI=[norminv(0.000005,xc_mean_(1),xc_std_(1)) norminv(0.999995,xc_mean_(1),xc_std_(1))];
axis([0.06 0.09 CI(1,1) CI(1,2)])
box on

figure(4)
scatter(x_v(:,1),x_v(:,3),'k.')
CI=[norminv(0.000005,xc_mean_(2),xc_std_(2)) norminv(0.999995,xc_mean_(2),xc_std_(2))];
axis([0.06 0.09 CI(1,1) CI(1,2)])
box on


