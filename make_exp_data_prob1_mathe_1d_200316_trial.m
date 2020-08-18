% 2020.02.22 Make Experiment Data
%%
clc
close all
clear
kkk=4;
for ttt=1:1000
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
dpt_n=[50 1;]; % dimension and number of experiment results
xc_mean_=15;xc_std_=1.3;

for nt=kkk%:4
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

folder_path=pwd;
cd(folder_path)

opt_folder=[pwd,'\opt_data_mathe_linear_200225'];
mkdir(opt_folder)

folder_name_type={'\exp_data_case1_normal_krig_200222';
    '\exp_data_case2_lognormal_krig_200222';
    '\exp_data_case3_gumbel_krig_200222';
    '\exp_data_case4_extreme_krig_200222';
    };
dist_type={'Normal','Lognormal','Extreme Value','Extreme Value'};
dist_name_type={'case1_normal';'case2_lognormal';'case3_gumbel';'case4_extreme';'case5_uniform';};
type_num=[1 2 3 6 5];
dpt_n=[50 1;]; % dimension and number of experiment results

%% Test - Input variable
% 시험 데이터 생성, disk, 190414
for nt=kkk%1:4
    % 저장
    folder_name='\exp_data_ex1_linear_200316';
    file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
        sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(2),dpt_n(1),1)];
    load(file_exp_name)
    x_sample(:,nt)=reshape(xc,size(xc,1)*size(xc,3),1);
end

for ph=kkk%1:4
    phat1 = mle(x_sample(:,ph),'distribution',dist_type{ph});
    phat2=phat1;
    if ph==2
        [M,N]=lognstat(phat1(1),phat1(2));
        phat2=[M,sqrt(N)];
    elseif ph==3
        phat1 = mle(-x_sample(:,ph),'distribution',dist_type{ph});
%         [M,N]=evstat(phat1(1),phat1(2));
        phat2=[-phat1(1)+0.5772*phat1(2),sqrt((-phat1(2))^2*pi^2/6)];
    elseif ph==4
        [M,N]=evstat(phat1(1),phat1(2));
        phat2=[M,sqrt(N)];
    elseif ph==5
        [M,N]=unifstat(phat1(1),phat1(2));
        phat2=[M,sqrt(N)];
    end
    parm_true2(ph,:)=phat2;
end


disp(sprintf('Trial%d',ttt))
if parm_true2(nt,1)>=14.80 && parm_true2(nt,1)<=15.20 && parm_true2(nt,2)>=1.285 && parm_true2(nt,2)<=1.315
    break
end
end