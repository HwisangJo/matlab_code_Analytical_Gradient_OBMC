%% Likelihood Function Metric
function [ff]=MLE_kerkrig_Ndim_analy_mathe_no_grad_corr_parfor(X_,file_exp_name,file_s_name,type_num,perf_func)
% ex) type_num=[1 3 5 2]; % for 4 unknown variables
% ex) X_=[mu1 sigma1 ... muN sigmaM]; % for N unknown variables
% 1: Normal, 2: Lognormal, 3: Gumbel, 5: Uniform, 6: Extreme
% perf_func=@(x,n_unknown) model_krig_10bar_ver02(x,n_unknown);
load(file_exp_name)
load(file_s_name)

dist_type={'Normal','Lognormal','Extreme Value',[],'Uniform','Extreme Value'};
dim=size(X_,2);
for h=1:dim/2
    if type_num(h)==5
        ab(1:2)=convert_param(X_((h-1)*2+1),X_(h*2),5);
        X((h-1)*2+1)=ab(1)-ab(2)*sqrt(3);
        X(h*2)=ab(1)+ab(2)*sqrt(3);
    else
        X((h-1)*2+1:h*2)=convert_param(X_((h-1)*2+1),X_(h*2),type_num(h));
    end
end
%% Calculate objective f
% MCS of input variable
n=10^6;
% n=10^5;
sim=zeros(n,1);
rng(s);
for h=1:dim/2
    if type_num(h)==3
        var_all(:,h)=X((h-1)*2+1)-(random(dist_type{type_num(h)},X((h-1)*2+1),X(h*2), [n, 1])-X((h-1)*2+1)); % Unknown Variables
    else
        var_all(:,h)=random(dist_type{type_num(h)},X((h-1)*2+1),X(h*2), [n, 1]); % Unknown Variables
    end
end

% % Performance Function
% perf_func=@(x) exp_sum(x);
% grad_func=@(x) exp_sum_grad(x);

% Uncertainty propagation(KDE) & Likelihood Function Metric
if nargout==1
for i=1:1:size(v_exp_data,2) %y
%     Xd=xd(:,i);
    Xd=[];
%     parfor j=1%:size(xd,1)
    for j=1%:size(xd,1)
%         y=perf_func([repmat(Xd(j,:),n,1), var_all],i);
        if size(v_exp_data,2)>1
        y=perf_func([var_all],i);
        elseif size(v_exp_data,2)==1
        y=perf_func([var_all]);
        end
        simm(:,j,i)=y;
    end
end
elseif nargout>1
for i=1:1:size(v_exp_data,2) %y
%     Xd=xd(:,i);
    Xd=[];
%     parfor j=1%:size(xd,1)
    for j=1%:size(xd,1)
%         [y,grad_y]=perf_func([repmat(Xd(j,:),n,1), var_all],i);
        if size(v_exp_data,2)>1
        [y,grad_y]=perf_func([var_all],i);
        elseif size(v_exp_data,2)==1
        [y,grad_y]=perf_func([var_all]);
        end
        simm(:,j,i)=y;
        Grad_y{j,i}=grad_y;
    end
end
end

d=size(v_exp_data,2);
m=(4/((d+2)*n))^(1/(d+4));
for j=1%:size(xd,1)
    for i=1:1:size(v_exp_data,2) %y
        mean_y(j,i)=mean(simm(:,j,i));
        sigma_y(j,i)=sqrt(sum((simm(:,j,i)-mean_y(j,i)).^2/(n-1)));
%         bw(j,i)=(4/(3*n))^(1/5)*sigma_y(j,i);
        bw(j,i)=m*sigma_y(j,i);
    end
    
    % KDE
    for ik=1:size(v_exp_data,3)
        for i=1:1:size(v_exp_data,2) %y
            pd_1(:,i)=normpdf(v_exp_data(j,i,ik),simm(:,j,i),bw(j,i));
        end
        pd_mv(ik,:)=sum(prod(pd_1,2),1)/n;
    end
    pd_(:,j)=pd_mv;
end

ff=sum(sum(-log(pd_)));
ff=ff/numel(pd_);

if ff>10^308
    ff=10^308;
end
% disp('%%%%%%%%%%%%%%%%%%%% function evaluation %%%%%%%%%%%%%%%%%%%%')

