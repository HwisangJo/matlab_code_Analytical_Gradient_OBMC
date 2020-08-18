% 2020.02.22 Make Experiment Data
%%
clc
close all
clear

addpath(genpath([pwd,'\imm1460\dace']))
folder_name_type={'\exp_data_case1_normal_krig_200222';
    '\exp_data_case2_lognormal_krig_200222';
    '\exp_data_case3_gumbel_krig_200222';
    '\exp_data_case4_extreme_krig_200222';
    };
dist_type={'Normal','Lognormal','Extreme Value','Extreme Value'};
dist_name_type={'case1_normal';'case2_lognormal';'case3_gumbel';'case4_extreme';'case5_uniform';};
type_num=[1 2 3 6 5];
%% 가상 시험 데이터 생성, 2020.02.22, Linear
% dpt_n=[6 9]; % dimension and number of experiment results
% dpt_n=[6 10]; % dimension and number of experiment results
dpt_n=[50 1;]; % dimension and number of experiment results
xc_mean_=15;xc_std_=1.0;
% xc_mean_=11;xc_std_=0.5;

for nt=1:4
    folder_name='\exp_data_ex1_linear_200316';
mkdir([pwd,folder_name,folder_name_type{nt}])
nexp=1;
for trial=1:nexp
    for pt=1
        % Observation Sites (6*10=60)
        ne=dpt_n(pt,2);tr=dpt_n(pt,1); % dimension and number of experiment results
        ab=convert_param(xc_mean_,xc_std_,type_num(nt));
        xc_mean=ab(1);xc_std=ab(2);
        if nt==3
            xc(:,1,:)=xc_mean-(random(dist_type{nt},xc_mean,xc_std, [ne,tr])-xc_mean); % Unknown Variables
        else
            xc(:,1,:)=random(dist_type{nt},xc_mean,xc_std, [ne,tr]); % Unknown Variables
        end
        
        % Model
        y_exp=[];x_exp=[];x_v=[];
        for i=1:tr
%             for k=1:ne
                v_exp_data(:,:,i)=model_ex_1d_linear([xc(:,:,i)]);
%                 v_exp_data(:,:,i)=model_krig([xd xc(:,:,i)]);
%                 v_exp_data(:,:,i)=branin([xd xc(:,:,i)]);
%             end
            x_v=[x_v;[xc(:,:,i)]];
            y_exp=[y_exp;v_exp_data(:,:,i)];
%             x_exp=[x_exp;[xd]];
        end
        
        % Sample mena, variance
        mean_samp=mean(x_v);
        sd_samp=std(x_v);
        var_samp=var(x_v);
        
        % 저장
        file_name=[pwd,folder_name,folder_name_type{nt},'\', ...
            sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(pt,2),dpt_n(pt,1),trial)];
        save(file_name,'v_exp_data','xc','ne','tr','y_exp','x_exp','mean_samp','sd_samp','var_samp','x_v') %
        clear v_exp_data xc ne tr y_exp x_exp x_v
    end
end
end

figure
% scatter(x_v,reshape(v_exp_data(1,1,:),30,1))
histogram(v_exp_data(1,1,:),8)

%% 가상 시험 데이터 생성, 2020.02.22, Non-linear
% dpt_n=[6 9]; % dimension and number of experiment results
% dpt_n=[6 10]; % dimension and number of experiment results
dpt_n=[50 1;]; % dimension and number of experiment results
xc_mean_=15;xc_std_=1.3;
% xc_mean_=11;xc_std_=0.5;

for nt=1:4
    folder_name='\exp_data_ex1_nonlinear_200316';
    folder_name_linear='\exp_data_ex1_linear_200316';
mkdir([pwd,folder_name,folder_name_type{nt}])
nexp=1;
for trial=1:nexp
    for pt=1
        file_name_linear=[pwd,folder_name_linear,folder_name_type{nt},'\', ...
            sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(pt,2),dpt_n(pt,1),trial)];
        load(file_name_linear)
%         % Observation Sites (6*10=60)
%         ne=dpt_n(pt,2);tr=dpt_n(pt,1); % dimension and number of experiment results
%         ab=convert_param(xc_mean_,xc_std_,type_num(nt));
%         xc_mean=ab(1);xc_std=ab(2);
%         if nt==3
%             xc(:,1,:)=xc_mean-(random(dist_type{nt},xc_mean,xc_std, [ne,tr])-xc_mean); % Unknown Variables
%         else
%             xc(:,1,:)=random(dist_type{nt},xc_mean,xc_std, [ne,tr]); % Unknown Variables
%         end
        
        % Model
        y_exp=[];x_exp=[];x_v=[];
        for i=1:tr
%             for k=1:ne
                v_exp_data(:,:,i)=model_ex_1d_nonlinear([xc(:,:,i)]);
%                 v_exp_data(:,:,i)=model_krig([xd xc(:,:,i)]);
%                 v_exp_data(:,:,i)=branin([xd xc(:,:,i)]);
%             end
            x_v=[x_v;[xc(:,:,i)]];
            y_exp=[y_exp;v_exp_data(:,:,i)];
%             x_exp=[x_exp;[xd]];
        end
        
        % Sample mena, variance
        mean_samp=mean(x_v);
        sd_samp=std(x_v);
        var_samp=var(x_v);
        
        % 저장
        file_name=[pwd,folder_name,folder_name_type{nt},'\', ...
            sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(pt,2),dpt_n(pt,1),trial)];
        save(file_name,'v_exp_data','xc','ne','tr','y_exp','x_exp','mean_samp','sd_samp','var_samp','x_v') %
        clear v_exp_data xc ne tr y_exp x_exp x_v
    end
end
end

figure
% scatter(x_v,reshape(v_exp_data(1,1,:),30,1))
histogram(v_exp_data(1,1,:),8)
