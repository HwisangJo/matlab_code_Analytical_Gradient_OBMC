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
%% 가상 시험 데이터 생성, 2020.02.22, 범위 넓은 조건
dpt_n=[50 1;]; % dimension and number of experiment results
% xc_mean_(1)=8000;xc_std_(1)=0.1*xc_mean_(1);
% xc_mean_(1)=210e9;xc_std_(1)=0.05*xc_mean_(1);
% xc_mean_(1)=210e9;xc_std_(1)=0.20*xc_mean_(1);
xc_mean_(1)=210e9;xc_std_(1)=0.10*xc_mean_(1);

perf_func=@(x) model_krig_modal_bracket2(x,1);
for nt=3%:4
    folder_name='\exp_data_ex2_bracket_1d_200316';
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
            v_exp_data(:,:,i)=perf_func([xc(:,:,i)]);
            x_v=[x_v;[xc(:,:,i)]];
            y_exp=[y_exp;v_exp_data(:,:,i)];
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

%%
for nt=4%1:2%:2
    for ini=1
        file_name=[pwd,folder_name,folder_name_type{nt},'\', ...
            sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(1,2),dpt_n(1,1),1)];
  load(file_name)      
type_num_=type_num(nt);
dist_type_={'Normal','Lognormal','Extreme Value',[],'Uniform','Extreme Value'};

CI(1,1)=norminv(0.00005,parm_true2(nt,1),parm_true2(nt,2));
CI(1,2)=norminv(0.99995,parm_true2(nt,1),parm_true2(nt,2));
sample_n=100;
% x1 = linspace(CI(1,1), CI(1,2),sample_n);
x1 = linspace(0, CI(1,2),sample_n);
xx=[x1'];

ab1=convert_param(parm_true2(nt,1),parm_true2(nt,2),type_num(nt));
xc_mean1=ab1(1);xc_std1=ab1(2);
if nt==3
    xc_mcs(:,1,:)=xc_mean1-(random(dist_type_{type_num_},xc_mean1,xc_std1, [10^6,1])-xc_mean1); % Unknown Variables
else
    xc_mcs(:,1,:)=random(dist_type_{type_num_},xc_mean1,xc_std1, [10^6,1]); % Unknown Variables
end

for iii=1%:2
ypdf_mle(:,iii)=pdf(dist_type_{type_num(nt)},xx(:,iii),xc_mean1,xc_std1);
if type_num(nt)==3
    ypdf_mle(:,iii)=pdf(dist_type_{type_num(nt)},-xx(:,iii)+2*parm_true2(nt,(iii-1)*2+1),parm_true2(nt,(iii-1)*2+1),parm_true2(nt,(iii-1)*2+2));
end
end

nbin1=20;
nbin2=1000;
Fig=figure(2*nt-1+10);
pp1=histogram(x_sample(:,nt),nbin1,'Normalization','pdf');
hold on
pp3=histogram(xc_mcs(:,1,1),nbin2,'Normalization','pdf');
% pp3=histogram(xc_mcs(:,1,1),nbin2);
hold on
pp2=plot(x1,ypdf_mle(:,1),'r-.','LineWidth',2.0);
% xlim([CI(1,1), CI(1,2)])
xlim([0, CI(1,2)])
% ylim([0 6]*10^-5)
set(gca,'fontsize', 13);
ytickangle(90)
xlabel('Youngs Modulus (\theta_{1})','fontsize', 15)
ylabel('PDF','fontsize', 15)
legend([pp1 pp2],'Virtual Input Data','mle PDF','Location','best','fontsize', 12)
set(Fig,'pos',[100 541 660 303]);
    end
end
