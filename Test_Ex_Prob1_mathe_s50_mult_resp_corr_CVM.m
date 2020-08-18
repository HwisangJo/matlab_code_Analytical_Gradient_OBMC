%% MLE inverse problem : KD-estimation (CVM)
clc
close all
clear

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
dpt_n=[50 1]; % dimension and number of experiment results
%% Test - Input variable
% 시험 데이터 생성, disk, 190414
for nt=1:2
    for kk=1:2
        % 저장
        file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
            sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200229.mat'],dpt_n(2),dpt_n(1),1)];
        load(file_exp_name)
        x_sample(:,kk,nt)=reshape(xc(:,kk,:),size(xc,1)*size(xc,3),1);
    end
end

for ph=1:2
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
foldername4= [opt_folder,'\','opt_input_mle'];
save(foldername4,'parm_true2') % sim?

%%
foldername5=[opt_folder,'\','opt_hist_cvm'];
mkdir(foldername5)
for nt=2%1:2
    for ini=1:5
        for si=1:6
            foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_norm_trr_cvm.mat', nt,ini,si)];
            
            % Optimization-Based Statistical Model Calibration (Parameter Estimation)
            file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
                sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200229.mat'],dpt_n(1,2),dpt_n(1,1),1)];
            file_s_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',si)];
            save_sign=0;
            
            % Performance Function
            perf_func=@(x,y_num) model_ex_mult_resp(x,y_num);
            
            lbx=[10, 10*0.05, 10, 10*0.05]*0.5;
            ubx=[10, 10*0.05, 10, 10*0.05]*1.5;

            ratio_x0=[1.10 1.05 0.90 0.95;0.92 0.98 1.10 1.10; 0.95 0.95 1.10 1.05;1.03 1.07 0.98 0.92;1.04 0.95 0.92 0.90];
            x0=parm_true2(nt,:).*ratio_x0(ini,:);
            
            lbx_n=zeros(1,size(lbx,2));
            ubx_n=ones(1,size(ubx,2));
            x0_n=(x0-lbx)./(ubx-lbx);
            
            obj = @(x) MLE_kerkrig_Ndim_numer_cvm_corr_parfor(x,10^-15,file_exp_name,file_s_name,type_num(nt,:),perf_func); % Likelihood Function Metric(KDE)
            obj_n = @(x) perf_norm(x,obj,lbx,ubx);
            
            % Optimization
            tic
            opts = optimoptions('fmincon','Display','iter','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'StepTolerance',1e-15,'OptimalityTolerance',1e-15,'ConstraintTolerance',1e-15,'MaxFunctionEvaluations',1000);
            [xvar_n,fval_sqp,history] = runfmincon(obj_n,x0_n,lbx_n,ubx_n);
            xvar=xvar_n.*(ubx-lbx)+lbx;
            XXX=[xvar(1) xvar(2) xvar(3) xvar(4) fval_sqp]
            time_opt_sqp(si,:)=toc
            
            save(foldername3,'xvar','fval_sqp','time_opt_sqp','lbx','ubx','history','x0') % sim?
            
            load([pwd,'\opt_data_mathe_mult_resp_corr50_200316\opt_input_mle.mat'])
            x_un=history.x.*(ubx-lbx)+lbx;
            
            Fig=figure(nt);
            set(Fig,'pos',[100 541 2260 743]);
            ax=subplot(1,2,1);
            scatter(x_un(1,1),x_un(1,2),'ro')
            hold on
            scatter(x_un(end,1),x_un(end,2),'rx')
            hold on
            plot(x_un(:,1),x_un(:,2),'k.-')
            box on
            hold on
            scatter(parm_true2(nt,1),parm_true2(nt,2),'bd')
            xlabel('Mean of \theta_{1}')
            ylabel('STD of \theta_{1}')
            title('History of TRR - \theta_{1}')
            xlim([lbx(1) ubx(1)])
            ylim([lbx(2) ubx(2)])
            
            ax=subplot(1,2,2);
            scatter(x_un(1,3),x_un(1,4),'ro')
            hold on
            scatter(x_un(end,3),x_un(end,4),'rx')
            hold on
            plot(x_un(:,3),x_un(:,4),'k.-')
            box on
            hold on
            scatter(parm_true2(nt,3),parm_true2(nt,4),'bd')
            xlabel('Mean of \theta_{2}')
            ylabel('STD of \theta_{2}')
            title('History of TRR - \theta_{2}')
            xlim([lbx(3) ubx(3)])
            ylim([lbx(4) ubx(4)])
            saveas(Fig,[foldername5,'\',sprintf('optimum_history_case%d_initial%d_seed%d',nt,ini,si),'.jpg'])
            close all
        end
    end
end
