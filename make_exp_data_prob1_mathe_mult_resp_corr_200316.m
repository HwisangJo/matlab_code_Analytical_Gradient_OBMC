% 2020.02.22 Make Experiment Data
%% design variable case1
clc
close all
clear

addpath(genpath([pwd,'\imm1460\dace']))
folder_name_type={'\exp_data_case1_normal_lognormal_200229';
    '\exp_data_case2_gumbel_extreme_200229';
    };
dist_type={'Normal','Lognormal';'Extreme Value','Extreme Value'};
dist_name_type={'case1_normal_lognormal';'case2_gumbel_extreme';};
type_num=[1 2; 3 6];

%% 가상 시험 데이터 생성 (Branin & Model sub - 15 design sites, same seeds with branin)
% dpt_n=[10 15;]; % dimension and number of experiment results
dpt_n=[50 1;]; % dimension and number of experiment results
xc_mean_=11;xc_std_=0.5;
xc_mean_(1)=10;xc_std_(1)=0.5;
xc_mean_(2)=10;xc_std_(2)=0.5;

for nt=1:2
folder_name='\exp_data_ex1_exp_multi_200316';
mkdir([pwd,folder_name,folder_name_type{nt}])
nexp=50;
for ei=1:nexp
    for pt=1
        file_s_exp_name=[pwd,'\s_save_exp_200228','\',sprintf('s_save_exp_trial%d_200228.mat',ei)];
        load(file_s_exp_name)
        rng(s);
        
        % Observation Sites (6*10=60)
        ne=dpt_n(pt,2);tr=dpt_n(pt,1); % dimension and number of experiment results
        for ui=1:2
        ab(ui,:)=convert_param(xc_mean_(ui),xc_std_(ui),type_num(nt,ui));
        xc_mean(ui)=ab(ui,1);xc_std(ui)=ab(ui,2);
        if type_num(nt,ui)==3
            xc(:,ui,:)=xc_mean(ui)-(random(dist_type{nt,ui},xc_mean(ui),xc_std(ui), [ne,tr])-xc_mean(ui)); % Unknown Variables
        else
            xc(:,ui,:)=random(dist_type{nt,ui},xc_mean(ui),xc_std(ui), [ne,tr]); % Unknown Variables
        end
        end
        
        % Model
        y_exp=[];x_exp=[];x_v=[];
        for ii=1:2
        for i=1:tr
%             v_exp_data(:,:,i)=branin([xd xc(:,:,i)]);
%             v_exp_data(:,:,i)=model_ex_sub([xd xc(:,:,i)]);
            v_exp_data(:,ii,i)=model_ex_mult_resp([xc(:,:,i)],ii);
            x_v=[x_v;[xc(:,:,i)]];
            y_exp=[y_exp;v_exp_data(:,ii,i)];
        end
        end
        
        % Sample mena, variance
        mean_samp=mean(x_v,1);
        sd_samp=std(x_v,1);
        var_samp=var(x_v,1);
        
        % 저장
        file_name=[pwd,folder_name,folder_name_type{nt},'\', ...
            sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200229.mat'],dpt_n(pt,2),dpt_n(pt,1),ei)];
        save(file_name,'v_exp_data','xc','ne','tr','y_exp','mean_samp','sd_samp','var_samp','x_v') %
        clear v_exp_data xc ne tr y_exp x_v
    end
end
end

scatter(v_exp_data(1,1,1:dpt_n(1)),v_exp_data(1,2,1:dpt_n(1)),'k.')

figure(2)
scatter(repmat(xd,tr,1),y_exp(:,1),'k.')
% axis([0 40 -20 190])
axis([0 10 30 140])
box on

figure(3)
scatter(x_v(:,1),x_v(:,2),'k.')
% axis([0 40 0 30])
axis([0 10 5 15])
box on
CI=[norminv(0.0005,xc_mean_,xc_std_) norminv(0.9995,xc_mean_,xc_std_)];

