% 2020.02.22 Make Experiment Data
%%
clc
close all
clear
for ttt=1:1000
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
xc_mean_(1)=210e9;xc_std_(1)=0.1*xc_mean_(1);
xc_mean_(2)=8000;xc_std_(2)=0.1*xc_mean_(2);

perf_func1=@(x) model_krig_modal_bracket(x,2);
perf_func2=@(x) model_krig_modal_bracket2(x,1);
for nt=2%:2
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
        save(file_name,'v_exp_data','xc','ne','tr','y_exp','mean_samp','sd_samp','var_samp','x_v') %
        clear v_exp_data xc ne tr y_exp x_v
    end
end
end

folder_path=pwd;
cd(folder_path)

% opt_folder=[pwd,'\opt_data_ex2_krig_200309'];
% mkdir(opt_folder)

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
for nt=2%1:2
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

disp(sprintf('Trial%d',ttt))
if parm_true2(nt,1)>=2.00*10^11 && parm_true2(nt,1)<=2.20*10^11 && parm_true2(nt,2)>=2.05*10^10 && parm_true2(nt,2)<=2.15*10^10 && ...
        parm_true2(nt,3)>=7.90*10^3 && parm_true2(nt,3)<=8.10*10^3 && parm_true2(nt,4)>=7.95*10^2 && parm_true2(nt,4)<=8.05*10^2
    break
end
end


