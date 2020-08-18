% 2020.02.22 Make Experiment Data
%% design variable case1
clc
close all
clear

for ttt=1:1000
addpath(genpath([pwd,'\imm1460\dace']))
folder_name_type={'\exp_data_case1_normal_lognormal_200229';
    '\exp_data_case2_gumbel_extreme_200229';
    };
dist_type={'Normal','Lognormal';'Extreme Value','Extreme Value'};
dist_name_type={'case1_normal_lognormal';'case2_gumbel_extreme';};
type_num=[1 2; 3 6];

%% 가상 시험 데이터 생성 (Branin & Model sub - 15 design sites, same seeds with branin)
dpt_n=[50 1;]; % dimension and number of experiment results
xc_mean_(1)=10;xc_std_(1)=0.5;
xc_mean_(2)=10;xc_std_(2)=0.5;

for nt=2%:2
folder_name='\exp_data_ex1_exp_multi_200316';
mkdir([pwd,folder_name,folder_name_type{nt}])
nexp=50;
for ei=1:nexp
    for pt=1
%         file_s_exp_name=[pwd,'\s_save_exp_200228','\',sprintf('s_save_exp_trial%d_200228.mat',ei)];
%         load(file_s_exp_name)
%         rng(s);
        
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


% folder_path='D:\Research\Research_for_paper1\Test_Ex_Prob1';
folder_path=pwd;
cd(folder_path)

opt_folder=[pwd,'\opt_data_mathe_mult_resp_corr50_200316'];
mkdir(opt_folder)

folder_name_type={'\exp_data_case1_normal_lognormal_200229';
    '\exp_data_case2_gumbel_extreme_200229';
    };
dist_type={'Normal','Lognormal';'Extreme Value','Extreme Value'};
dist_name_type={'case1_normal_lognormal';'case2_gumbel_extreme';};
type_num=[1 2; 3 6];
folder_name='\exp_data_ex1_exp_multi_200316';
% dpt_n=[10 15;]; % dimension and number of experiment results
dpt_n=[50 1]; % dimension and number of experiment results
%% Test - Input variable
% 시험 데이터 생성, disk, 190414
for nt=2%:2
    for kk=1:2
    % 저장
    file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
    sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200229.mat'],dpt_n(2),dpt_n(1),1)];
    load(file_exp_name)
    x_sample(:,kk,nt)=reshape(xc(:,kk,:),size(xc,1)*size(xc,3),1);
    end
end

for ph=2%:2
    for kk=1:2
    phat1 = mle(x_sample(:,kk,ph),'distribution',dist_type{ph,kk});
    phat2=phat1;
    if type_num(ph,kk)==2
        [M,N]=lognstat(phat1(1),phat1(2));
        phat2=[M,sqrt(N)];
    elseif type_num(ph,kk)==3
        phat1 = mle(-x_sample(:,kk,ph),'distribution',dist_type{ph});
%         [M,N]=evstat(phat1(1),phat1(2));
        phat2=[-phat1(1)+0.5772*phat1(2),sqrt((-phat1(2))^2*pi^2/6)];
    elseif type_num(ph,kk)==6
        [M,N]=evstat(phat1(1),phat1(2));
        phat2=[M,sqrt(N)];
    elseif type_num(ph,kk)==5
        [M,N]=unifstat(phat1(1),phat1(2));
        phat2=[M,sqrt(N)];
    end
    parm_true2(ph,2*(kk-1)+[1:2])=phat2;
    end
end

disp(sprintf('Trial%d',ttt))
if parm_true2(nt,1)>=9.70 && parm_true2(nt,1)<=10.30 && parm_true2(nt,2)>=0.470 && parm_true2(nt,2)<=0.530 && ...
        parm_true2(nt,3)>=9.70 && parm_true2(nt,3)<=10.30 && parm_true2(nt,4)>=0.470 && parm_true2(nt,4)<=0.530
    break
end
end
