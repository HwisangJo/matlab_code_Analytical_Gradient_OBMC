clc
close all
clear

%% Error Calculation
sss=0;
opt_folder=[pwd,'\opt_data_ex2_krig_200309'];
% dist_name_type={'case1_normal_gumbel';'case2_lognormal_extreme'};
dist_name_type={'case1_normal_lognormal';'case2_gumbel_extreme'};
foldername4= [opt_folder,'\','opt_input_mle'];
load(foldername4)
for nt=1%:2
    for ini=1:5
        for i=1:6
    foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_norm_trr_analy.mat', nt,ini,i)];
    load(foldername3)
    err(nt,:,ini,i)=abs((parm_true2(nt,:)-xvar))./parm_true2(nt,:)*100;
    parm_xvar(nt,:,ini,i)=xvar;
%     n_rmse(nt,:,i)=sqrt(sum((parm_true2(nt,:)-parm_xvar(nt,:,i)).^2,3)/size(parm_xvar(nt,:,i),3))./parm_true2(nt,:)*100;
    if fval_sqp==10^308
        sss=sss+1;
        save_ind(sss,:)=[nt ini i];
        disp('initial값이 무한으로')
    end
        end
    end
    r_rmse_analy(nt,:)=sqrt(sum(sum(((parm_true2(nt,:)-parm_xvar(nt,:,:,:))).^2,4),3)/(size(parm_xvar,3)*size(parm_xvar,4)))./parm_true2(nt,:)*100;
%     r_rmse(nt,:)=sqrt(sum(sum(((parm_true2(nt,:)-parm_xvar(nt,:,[1 3],:))).^2,4),3)/12)./parm_true2(nt,:)*100;
end
parm_xvar_table_analy=sum(sum(parm_xvar,4),3)/(size(parm_xvar,3)*size(parm_xvar,4));

sss=0;
for nt=1:2
    for ini=1:5
        for i=1:6
    foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_norm_trr_fdm.mat', nt,ini,i)];
    load(foldername3)
    err(nt,:,ini,i)=abs((parm_true2(nt,:)-xvar))./parm_true2(nt,:)*100;
    parm_xvar(nt,:,ini,i)=xvar;
%     n_rmse(nt,:,i)=sqrt(sum((parm_true2(nt,:)-parm_xvar(nt,:,i)).^2,3)/size(parm_xvar(nt,:,i),3))./parm_true2(nt,:)*100;
    if fval_sqp==10^308
        sss=sss+1;
        save_ind(sss,:)=[nt ini i];
        disp('initial값이 무한으로')
    end
        end
    end
    r_rmse_fdm(nt,:)=sqrt(sum(sum(((parm_true2(nt,:)-parm_xvar(nt,:,:,:))).^2,4),3)/(size(parm_xvar,3)*size(parm_xvar,4)))./parm_true2(nt,:)*100;
end
parm_xvar_table_fdm=mean(mean(parm_xvar,4),3);

sss=0;
for nt=1:2
    for ini=1:5
        for i=1:6
    foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_norm_trr_cvm_10_4.mat', nt,ini,i)];
    load(foldername3)
    err(nt,:,ini,i)=abs((parm_true2(nt,:)-xvar))./parm_true2(nt,:)*100;
    parm_xvar(nt,:,ini,i)=xvar;
%     n_rmse(nt,:,i)=sqrt(sum((parm_true2(nt,:)-parm_xvar(nt,:,i)).^2,3)/size(parm_xvar(nt,:,i),3))./parm_true2(nt,:)*100;
    if fval_sqp==10^308
        sss=sss+1;
        save_ind(sss,:)=[nt ini i];
        disp('initial값이 무한으로')
    end
        end
    end
    r_rmse_cvm_10_4(nt,:)=sqrt(sum(sum(((parm_true2(nt,:)-parm_xvar(nt,:,:,:))).^2,4),3)/(size(parm_xvar,3)*size(parm_xvar,4)))./parm_true2(nt,:)*100;
end
parm_xvar_table_cvm_10_4=mean(mean(parm_xvar,4),3);

sss=0;
for nt=1:2
    for ini=1:5
        for i=1:6
    foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_norm_trr_analy_10_4.mat', nt,ini,i)];
    load(foldername3)
    err(nt,:,ini,i)=abs((parm_true2(nt,:)-xvar))./parm_true2(nt,:)*100;
    parm_xvar(nt,:,ini,i)=xvar;
%     n_rmse(nt,:,i)=sqrt(sum((parm_true2(nt,:)-parm_xvar(nt,:,i)).^2,3)/size(parm_xvar(nt,:,i),3))./parm_true2(nt,:)*100;
    if fval_sqp==10^308
        sss=sss+1;
        save_ind(sss,:)=[nt ini i];
        disp('initial값이 무한으로')
    end
        end
    end
    r_rmse_analy_10_4(nt,:)=sqrt(sum(sum(((parm_true2(nt,:)-parm_xvar(nt,:,:,:))).^2,4),3)/(size(parm_xvar,3)*size(parm_xvar,4)))./parm_true2(nt,:)*100;
end
parm_xvar_table_analy_10_4=mean(mean(parm_xvar,4),3);