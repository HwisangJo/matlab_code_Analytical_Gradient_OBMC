clc
close all
clear

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
for nt=1:4
    % 저장
    folder_name='\exp_data_ex1_linear_200316';
    file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
        sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(2),dpt_n(1),1)];
    load(file_exp_name)
    x_sample(:,nt)=reshape(xc,size(xc,1)*size(xc,3),1);
end

for ph=1:4
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
foldername4= [opt_folder,'\','opt_input_mle'];
save(foldername4,'parm_true2') % sim?

%% Estimated Input 1D
for nt=1%:4
    % 저장
    file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
        sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(2),dpt_n(1),1)];
    load(file_exp_name)
    x_sample(:,nt)=reshape(xc,size(xc,1)*size(xc,3),1);
end

foldername4= [opt_folder,'\','opt_input_mle'];
load(foldername4)
pt=1;si=1;
for nt=1%:4
    for ini=1
        foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_analy.mat', nt,ini,si)];
        load(foldername3)
        
        type_num_=type_num(nt);
        dist_type_={'Normal','Lognormal','Extreme Value',[],'Uniform','Extreme Value'};
        X_=xvar;
        dim=size(X_,2);
        for h=1:dim/2
            if type_num_(h)==5
                ab(1:2)=convert_param(X_((h-1)*2+1),X_(h*2),5);
                X((h-1)*2+1)=ab(1)-ab(2)*sqrt(3);
                X(h*2)=ab(1)+ab(2)*sqrt(3);
            else
                X((h-1)*2+1:h*2)=convert_param(X_((h-1)*2+1),X_(h*2),type_num_(h));
            end
        end
        
        ratio_x0=[1.08 1.07 0.90 0.95;0.92 0.98 1.10 1.10; 0.95 0.95 1.10 1.05;1.03 1.07 0.98 0.92;1.04 0.95 0.92 0.90];
        X0=parm_true2(nt,:).*ratio_x0(ini,1:2);
        for h=1:dim/2
            if type_num_(h)==5
                ab(1:2)=convert_param(X0((h-1)*2+1),X0(h*2),5);
                X0_((h-1)*2+1)=ab(1)-ab(2)*sqrt(3);
                X0_(h*2)=ab(1)+ab(2)*sqrt(3);
            else
                X0_((h-1)*2+1:h*2)=convert_param(X0((h-1)*2+1),X0(h*2),type_num_(h));
            end
        end
        
        CI(1,1)=norminv(0.00005,xvar(1),xvar(2));CI(1,2)=norminv(0.99995,xvar(1),xvar(2));
        CI_(1,1)=norminv(0.00005,X0(1),X0(2));CI_(1,2)=norminv(0.99995,X0(1),X0(2));
        sample_n=100;
        x1 = linspace(min(CI(1,1),CI_(1,1)), max(CI(1,2),CI_(1,2)),sample_n);
        xx=[x1'];
        
        for iii=1%:2
            ypdf_ini(:,iii)=pdf(dist_type_{type_num(nt,iii)},xx(:,iii),X0_((iii-1)*2+1),X0_((iii-1)*2+2));
            ypdf_cal(:,iii)=pdf(dist_type_{type_num(nt,iii)},xx(:,iii),X((iii-1)*2+1),X((iii-1)*2+2));
            if type_num_(:,iii)==3
                ypdf_ini(:,iii)=pdf(dist_type_{type_num(nt,iii)},-xx(:,iii)+2*X0_((iii-1)*2+1),X0_((iii-1)*2+1),X0_((iii-1)*2+2));
                ypdf_cal(:,iii)=pdf(dist_type_{type_num(nt,iii)},-xx(:,iii)+2*X((iii-1)*2+1),X((iii-1)*2+1),X((iii-1)*2+2));
            end
        end
        
        nbin1=5;
        Fig=figure(2*nt-1);
        pp1=histogram(x_sample(:,nt),nbin1,'Normalization','pdf','FaceColor',[0.3010 0.7450 0.9330]);
        hold on
        pp2=plot(x1,ypdf_ini(:,1),'r-.','LineWidth',2.0);
        hold on
        pp3=plot(x1,ypdf_cal(:,1),'k-','LineWidth',2.0);
        xlim([min(CI(1,1),CI_(1,1)), max(CI(1,2),CI_(1,2))])
        ylim([0 4]*10^-1)
        set(gca,'fontsize', 15);
        ytickangle(90)
        xlabel('\theta','fontsize', 17)
        ylabel('PDF','fontsize', 17)
        legend([pp1 pp2 pp3],'Virtual input data','Prior PDF','Estimated PDF','Location','best','fontsize', 12)
        set(Fig,'pos',[100 541 660 303]);
        set(gca,'fontname','times')  % Set it to times
    end
end

%% Output PDF Plot
foldername5=[opt_folder,'\','opt_hist_analy'];
mkdir(foldername5)
for nt=1%:4
    % poolobj=parpool(8);
    for ini=1%:5
        for si=1%:6
            foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_analy.mat', nt,ini,si)];
            
            % Optimization-Based Statistical Model Calibration (Parameter Estimation)
            folder_name='\exp_data_ex1_nonlinear_200316';
            file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
                sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200222.mat'],dpt_n(2),dpt_n(1),1)];
            file_s_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',si)];
            save_sign=0;
            
            % Performance Function
            perf_func=@(x) model_ex_1d_nonlinear(x);
            
            ratio_x0=[1.08 1.07 0.90 0.95;0.92 0.98 1.10 1.10; 0.95 0.95 1.10 1.05;1.03 1.07 0.98 0.92;1.04 0.95 0.92 0.90];
            x0=parm_true2(nt,:).*ratio_x0(ini,1:2);
            
            load(foldername3)
        end
    end
end

% min(v_exp_data,[],3)
% max(v_exp_data,[],3)

y_range1_=[4.8 6.2];
y_range1=[4.8 6.2];

load(file_s_name)
load(file_exp_name)

n=10^6;
sim=zeros(n,1);
rng(s);
dim=size(x0,2);
X=x0;
type_num_=type_num(1,:);
dist_type_={'Normal','Lognormal','Extreme Value',[],'Uniform','Extreme Value'};
for h=1:dim/2
    if type_num_(:,h)==3
        var_all(:,h)=X((h-1)*2+1)-(random(dist_type_{type_num_(:,h)},X((h-1)*2+1),X(h*2), [n, 1])-X((h-1)*2+1)); % Unknown Variables
    else
        var_all(:,h)=random(dist_type_{type_num_(:,h)},X((h-1)*2+1),X(h*2), [n, 1]); % Unknown Variables
    end
end
for i=1:1:size(v_exp_data,2) %y
    %     Xd=xd(:,i);
    Xd=[];
    %     parfor j=1%:size(xd,1)
    for j=1%:size(xd,1)
        %         y=perf_func([repmat(Xd(j,:),n,1), var_all],i);
        y=perf_func([var_all]);
        simm(:,j,i)=y;
    end
end

d=size(v_exp_data,2);
m=(4/((d+2)*n))^(1/(d+4));
j=1;
for i=1:1:size(v_exp_data,2) %y
    mean_y(j,i)=mean(simm(:,j,i));
    sigma_y(j,i)=sqrt(sum((simm(:,j,i)-mean_y(j,i)).^2/(n-1)));
    %         bw(j,i)=(4/(3*n))^(1/5)*sigma_y(j,i);
    bw(j,i)=m*sigma_y(j,i);
end

X_=linspace(y_range1_(1),y_range1_(2),500)';
for t=1:size(X_,1)
    for i=1:1:size(v_exp_data,2) %y
        pd_1(:,t)=normpdf(X_(t,i),simm(:,1,i),bw(1,i));
    end
end
pd__=(sum(prod(pd_1,3),1)/size(simm,1))';


n=10^6;
sim=zeros(n,1);
rng(s);
dim=size(x0,2);
X=xvar;
type_num_=type_num(1,:);
dist_type_={'Normal','Lognormal','Extreme Value',[],'Uniform','Extreme Value'};
for h=1:dim/2
    if type_num_(:,h)==3
        var_all(:,h)=X((h-1)*2+1)-(random(dist_type_{type_num_(:,h)},X((h-1)*2+1),X(h*2), [n, 1])-X((h-1)*2+1)); % Unknown Variables
    else
        var_all(:,h)=random(dist_type_{type_num_(:,h)},X((h-1)*2+1),X(h*2), [n, 1]); % Unknown Variables
    end
end
for i=1:1:size(v_exp_data,2) %y
    %     Xd=xd(:,i);
    Xd=[];
    %     parfor j=1%:size(xd,1)
    for j=1%:size(xd,1)
        y=perf_func([var_all]);
        simm(:,j,i)=y;
    end
end

d=size(v_exp_data,2);
m=(4/((d+2)*n))^(1/(d+4));
j=1;
for i=1:1:size(v_exp_data,2) %y
    mean_y(j,i)=mean(simm(:,j,i));
    sigma_y(j,i)=sqrt(sum((simm(:,j,i)-mean_y(j,i)).^2/(n-1)));
    %         bw(j,i)=(4/(3*n))^(1/5)*sigma_y(j,i);
    bw(j,i)=m*sigma_y(j,i);
end

X=linspace(y_range1(1),y_range1(2),500)';
for t=1:size(X,1)
    for i=1:1:size(v_exp_data,2) %y
        pd_2(:,t)=normpdf(X(t,i),simm(:,1,i),bw(1,i));
    end
end
pd_=(sum(prod(pd_2,3),1)/size(simm,1))';

nbin1=20;
Fig=figure;
set(Fig,'pos',[100 541 660 303]);
hold on
s3=histogram(v_exp_data(1,1,1:50),nbin1,'Normalization','pdf','FaceColor',[0.9290 0.6940 0.1250]);
hold on
s2=plot(X_,pd__,'r-.','LineWidth',2.0);
hold on
s1=plot(X,pd_,'k-','LineWidth',1.5);
box on
xlim(y_range1_)
% ylim([])
% grid on
legend([s2 s1 s3],{'Initial output PDF','Calibrated output PDF','Experimental data'})
set(gca,'fontname','times','fontsize', 15)  % Set it to times
xlabel('\it y','fontsize', 17)
ylabel('PDF','fontsize', 17)
axx = gca;
% axx.FontSize = 20;