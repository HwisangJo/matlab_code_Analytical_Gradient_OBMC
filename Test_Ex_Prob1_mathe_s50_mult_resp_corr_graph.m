%% MLE inverse problem : KD-estimation
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

%% Contour Plot
y_range1_=[20 45];
y_range2_=[9 18];
y_range1=[30 70];
y_range2=[7 13];

for nt=1%:2
    for ini=1%:5
        for si=1%:6
            foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_norm_trr_analy.mat', nt,ini,si)];
            
            % Optimization-Based Statistical Model Calibration (Parameter Estimation)
            file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
                sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200229.mat'],dpt_n(1,2),dpt_n(1,1),1)];
            file_s_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',si)];
            % file_s_exp_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',ei)];
            % load(file_s_exp_name)
            save_sign=0;
            
            % Performance Function
            perf_func=@(x,y_num) model_ex_mult_resp(x,y_num);
            
            lbx=[10, 10*0.05, 10, 10*0.05]*0.5;
            ubx=[10, 10*0.05, 10, 10*0.05]*1.5;
            ratio_x0=[1.10 1.05 0.90 0.95;0.92 0.98 1.10 1.10; 0.95 0.95 1.10 1.05;1.03 1.07 0.98 0.92;1.04 0.95 0.92 0.90];
            x0=parm_true2(nt,:).*ratio_x0(ini,:);
            
            load(foldername3) % sim?
        end
    end
end

load(file_s_name)
load(file_exp_name)

n=10^6;
sim=zeros(n,1);
rng(s);
dim=size(x0,2);
X=x0;
type_num_=type_num(2,:);
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
        y=perf_func([var_all],i);
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

[X_,Y_] = meshgrid(linspace(y_range1_(1),y_range1_(2),50),linspace(y_range2_(1),y_range2_(2),50));
for k=1:size(X_,1)
    grid_=[X_(:,k) Y_(:,k)];
    for t=1:size(X_,1)
        for i=1:1:size(v_exp_data,2) %y
            pd_1(:,t,i)=normpdf(grid_(t,i),simm(:,1,i),bw(1,i));
        end
    end
    pd__(k,:)=sum(prod(pd_1,3),1)/size(simm,1);
    disp(sprintf('%d th complete',k))
end

n=10^6;
sim=zeros(n,1);
rng(s);
dim=size(x0,2);
X=xvar;
type_num_=type_num(2,:);
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
        y=perf_func([var_all],i);
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

[X,Y] = meshgrid(linspace(y_range1(1),y_range1(2),50),linspace(y_range2(1),y_range2(2),50));
for k=1:size(X,1)
    grid_=[X(:,k) Y(:,k)];
    for t=1:size(X,1)
        for i=1:1:size(v_exp_data,2) %y
            pd_1(:,t,i)=normpdf(grid_(t,i),simm(:,1,i),bw(1,i));
        end
    end
    pd_(k,:)=sum(prod(pd_1,3),1)/size(simm,1);
    disp(sprintf('%d th complete',k))
end

% [f,xi]=ksdensity([simm(:,1,1) simm(:,1,2)],[X(:),Y(:)],'Bandwidth',bw);
% ff=reshape(f,50,50);

% scatter(simm(:,1,1),simm(:,1,2),'k.')
Fig=figure;
set(Fig,'pos',[700 541 710 593]);
hold on
s2=contour(X_,Y_,pd__',':','LineWidth',2.0);
hold on
s1=contour(X,Y,pd_','-','LineWidth',1.5);
% hold on
% s2=contour(X,Y,ff,':','LineWidth',2.0);
hold on
% scatter(v_exp_data(1,1,1:50),v_exp_data(1,2,1:50),'rx')
s3=scatter(v_exp_data(1,1,1:50),v_exp_data(1,2,1:50),100,'k.');
box on
axx = gca;
axx.FontSize = 17;
xlim([20 70])
ylim([7 18])
xlabel('\it y_{1}','fontsize',20)
ylabel('\it y_{2}','fontsize',20)
grid on
legend({'Initial output PDF','Calibrated output PDF','Experimental data'},'Location','best','fontsize',17)
set(gca,'fontname','times')  % Set it to times

%% Estimated Input 1D
for nt=1:2
    for kk=1:2
        % 저장
        file_exp_name=[pwd,folder_name,folder_name_type{nt},'\', ...
            sprintf(['exp_data_ex1_',dist_name_type{nt},'_n%d_%d_trial%d_200229.mat'],dpt_n(2),dpt_n(1),1)];
        load(file_exp_name)
        x_sample(:,kk,nt)=reshape(xc(:,kk,:),size(xc,1)*size(xc,3),1);
    end
end

foldername4= [opt_folder,'\','opt_input_mle'];
load(foldername4)
pt=1;si=1;
for nt=1%:2
    for ini=1
        foldername3= [opt_folder,'\','opt_',dist_name_type{nt},sprintf('_inform%d_%d_%d_%d.mat', nt,dpt_n(pt,2),si,ini)];
        load(foldername3)
        
        type_num_=type_num(nt,:);
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
        
        ratio_x0=[1.10 1.05 0.90 0.95;0.92 0.98 1.10 1.10; 0.95 0.95 1.10 1.05;1.03 1.07 0.98 0.92;1.04 0.95 0.92 0.90];
        X0=parm_true2(nt,:).*ratio_x0(ini,:);
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
        CI(2,1)=norminv(0.00005,xvar(3),xvar(4));CI(2,2)=norminv(0.99995,xvar(3),xvar(4));
        CI_(1,1)=norminv(0.00005,X0(1),X0(2));CI_(1,2)=norminv(0.99995,X0(1),X0(2));
        CI_(2,1)=norminv(0.00005,X0(3),X0(4));CI_(2,2)=norminv(0.99995,X0(3),X0(4));
        sample_n=100;
        x1 = linspace(min(CI(1,1),CI_(1,1)), max(CI(1,2),CI_(1,2)),sample_n);
        x2 = linspace(min(CI(2,1),CI_(2,1)), max(CI(2,2),CI_(2,2)),sample_n);
        xx=[x1' x2'];
        
        for iii=1:2
            ypdf_ini(:,iii)=pdf(dist_type_{type_num(nt,iii)},xx(:,iii),X0_((iii-1)*2+1),X0_((iii-1)*2+2));
            ypdf_cal(:,iii)=pdf(dist_type_{type_num(nt,iii)},xx(:,iii),X((iii-1)*2+1),X((iii-1)*2+2));
            if type_num_(:,iii)==3
                ypdf_ini(:,iii)=pdf(dist_type_{type_num(nt,iii)},-xx(:,iii)+2*X0_((iii-1)*2+1),X0_((iii-1)*2+1),X0_((iii-1)*2+2));
                ypdf_cal(:,iii)=pdf(dist_type_{type_num(nt,iii)},-xx(:,iii)+2*X((iii-1)*2+1),X((iii-1)*2+1),X((iii-1)*2+2));
            end
        end
        
        nbin1=8;
        Fig=figure(2*nt-1);
        pp1=histogram(x_sample(:,1,nt),nbin1,'Normalization','pdf');
        hold on
        pp2=plot(x1,ypdf_ini(:,1),'r-.','LineWidth',2.0);
        hold on
        pp3=plot(x1,ypdf_cal(:,1),'k-','LineWidth',2.0);
        xlim([min(CI(1,1),CI_(1,1)), max(CI(1,2),CI_(1,2))])
        % ylim([0 6]*10^-5)
        set(gca,'fontsize', 13);
        ytickangle(90)
        xlabel('\theta_{1}','fontsize', 15)
        ylabel('PDF','fontsize', 15)
        legend([pp1 pp2 pp3],'Virtual input data','Prior PDF','Estimated PDF','Location','best','fontsize', 12)
        set(Fig,'pos',[100 541 660 303]);
        set(gca,'fontname','times')  % Set it to times
        
        nbin2=8;
        % nbin2=nbin1;
        Fig=figure(2*nt);
        pp1=histogram(x_sample(:,2,nt),nbin2,'Normalization','pdf');
        hold on
        pp2=plot(x2,ypdf_ini(:,2),'r-.','LineWidth',2.0);
        hold on
        pp3=plot(x2,ypdf_cal(:,2),'k-','LineWidth',2.0);
        xlim([min(CI(2,1),CI_(2,1)), max(CI(2,2),CI_(2,2))])
        % ylim([0 7]*10^-4)
        set(gca,'fontsize', 13);
        ytickangle(90)
        xlabel('\theta_{2}','fontsize', 15)
        ylabel('PDF','fontsize', 15)
        % legend([pp1 pp2 pp3],'Virtual Input Data','Prior PDF','Calibrated PDF','Location','best','fontsize', 12)
        box on
        set(Fig,'pos',[100 541 660 303]);
        set(gca,'fontname','times')  % Set it to times
    end
end
