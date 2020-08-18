clc
close all
clear

%% Error Calculation
opt_folder=[pwd,'\opt_data_mathe_bracket_1d_200316'];
dist_name_type={'case1_normal';'case2_lognormal';'case3_gumbel';'case4_extreme';'case5_uniform';};
foldername4= [opt_folder,'\','opt_input_mle'];
load(foldername4)
for nt=1:4
    for ini=1:5
        for si=1:6
            foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_analy.mat', nt,ini,si)];
            load(foldername3)
            err(nt,:,si)=abs((parm_true2(nt,:)-xvar))./parm_true2(nt,:)*100;
            parm_xvar(nt,:,ini,si)=xvar;
        end
        n_rmse(nt,:)=sqrt(sum((parm_true2(nt,:)-parm_xvar(nt,:,:)).^2,3)/size(parm_xvar(nt,:,:),3))./parm_true2(nt,:)*100;
    end
%     r_rmse(nt,:)=sqrt(sum(sum(((parm_true2(nt,:)-parm_xvar(nt,:,:,:))).^2,4),3)/(size(parm_xvar,3)*size(parm_xvar,4)))./parm_true2(nt,:)*100;
    r_rmse_analy(nt,:)=sqrt(sum(sum(((parm_true2(nt,:)-parm_xvar(nt,:,:,:))).^2,4),3)/(size(parm_xvar,3)*size(parm_xvar,4)))./parm_true2(nt,:)*100;
end
parm_xvar_table_analy=sum(sum(parm_xvar,4),3)/(size(parm_xvar,3)*size(parm_xvar,4));

opt_folder=[pwd,'\opt_data_mathe_bracket_1d_200316'];
dist_name_type={'case1_normal';'case2_lognormal';'case3_gumbel';'case4_extreme';'case5_uniform';};
foldername4= [opt_folder,'\','opt_input_mle'];
load(foldername4)
for nt=1:4
    for ini=1:5
        for si=1:6
            foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_cvm.mat', nt,ini,si)];
            load(foldername3)
            err(nt,:,si)=abs((parm_true2(nt,:)-xvar))./parm_true2(nt,:)*100;
            parm_xvar(nt,:,ini,si)=xvar;
        end
        n_rmse(nt,:)=sqrt(sum((parm_true2(nt,:)-parm_xvar(nt,:,:)).^2,3)/size(parm_xvar(nt,:,:),3))./parm_true2(nt,:)*100;
    end
    r_rmse_cvm(nt,:)=sqrt(sum(sum((parm_true2(nt,:)-parm_xvar(nt,:,:,:)).^2,4),3)/(size(parm_xvar,3)*size(parm_xvar,4)))./parm_true2(nt,:)*100;
end
parm_xvar_table_cvm=sum(sum(parm_xvar,4),3)/(size(parm_xvar,3)*size(parm_xvar,4));

opt_folder=[pwd,'\opt_data_mathe_bracket_1d_200316'];
dist_name_type={'case1_normal';'case2_lognormal';'case3_gumbel';'case4_extreme';'case5_uniform';};
foldername4= [opt_folder,'\','opt_input_mle'];
load(foldername4)
for nt=1:4
    for ini=1:5
        for si=1:6
            foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_fdm.mat', nt,ini,si)];
            load(foldername3)
            err(nt,:,si)=abs((parm_true2(nt,:)-xvar))./parm_true2(nt,:)*100;
            parm_xvar(nt,:,ini,si)=xvar;
        end
        n_rmse(nt,:)=sqrt(sum((parm_true2(nt,:)-parm_xvar(nt,:,:)).^2,3)/size(parm_xvar(nt,:,:),3))./parm_true2(nt,:)*100;
    end
    r_rmse_fdm(nt,:)=sqrt(sum(sum((parm_true2(nt,:)-parm_xvar(nt,:,:,:)).^2,4),3)/(size(parm_xvar,3)*size(parm_xvar,4)))./parm_true2(nt,:)*100;
end
parm_xvar_table_fdm=sum(sum(parm_xvar,4),3)/(size(parm_xvar,3)*size(parm_xvar,4));

opt_folder=[pwd,'\opt_data_mathe_bracket_1d_200316'];
dist_name_type={'case1_normal';'case2_lognormal';'case3_gumbel';'case4_extreme';'case5_uniform';};
foldername4= [opt_folder,'\','opt_input_mle'];
load(foldername4)
for nt=1%:4
    for ini=1:5
        for si=1:6
            foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_youn_linear_normal.mat', nt,ini,si)];
            load(foldername3)
            err(nt,:,si)=abs((parm_true2(nt,:)-xvar))./parm_true2(nt,:)*100;
            parm_xvar(nt,:,ini,si)=xvar;
        end
        n_rmse(nt,:)=sqrt(sum((parm_true2(nt,:)-parm_xvar(nt,:,:)).^2,3)/size(parm_xvar(nt,:,:),3))./parm_true2(nt,:)*100;
    end
    r_rmse_linear(nt,:)=sqrt(sum(sum((parm_true2(nt,:)-parm_xvar(nt,:,:,:)).^2,4),3)/(size(parm_xvar,3)*size(parm_xvar,4)))./parm_true2(nt,:)*100;
end
parm_xvar_table_linear=sum(sum(parm_xvar,4),3)/(size(parm_xvar,3)*size(parm_xvar,4));