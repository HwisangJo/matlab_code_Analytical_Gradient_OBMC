%% Likelihood Function Metric
function [ff, gg_mult,time_analy]=MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(X_,file_exp_name,file_s_name,type_num,perf_func)
% function [ff, gg_mult]=MLE_kerkrig_Ndim_analy_mathe_multi_resp_corr_parfor(X_,file_exp_name,file_s_name,type_num,perf_func)

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
if nargout > 1 % gradient required
    tic
%% gradient required
N=size(xc,3);
for t=1%:size(xd,1) % 1 for paper
for i=1:1:size(v_exp_data,2)
    df=2*(simm(:,t,i)-reshape(v_exp_data(t,i,:),1,size(v_exp_data,3))); %l
    F=(0.5*df).^2;
    G=2*bw(t,i)^2;
    Af(:,:,i)=exp(-F/G);
    Bf(:,:,i)=F/G^2;
    Cf(:,:,i)=df/G;

    sim_dev(:,i,t)=simm(:,t,i)-mean_y(t,i);
    df1_2(:,i,t)=(N/((n-1)*sigma_y(t,i).^2))*sim_dev(:,i,t); % 2
    dg(:,i,t)=4*sim_dev(:,i,t)/(n-1)*m^2;
end
Af_t=prod(Af,3);
AAf=sum(Af_t,1);
for i=1:1:size(v_exp_data,2)
    df1_1_(:,:,i)=Af_t.*Bf(:,:,i);
    df1_1=sum(df1_1_,1);
    ddf1_11(:,:,i)=df1_1(:,:,i).*dg(:,i,t);
    df1_12(:,:,i)=-Af_t.*Cf(:,:,i);
    ddf1(:,i,t)=sum(-(ddf1_11(:,:,i)+df1_12(:,:,i))./AAf,2)+df1_2(:,i,t);
end
end

for i=1:1:size(v_exp_data,2)
    for h=1:dim/2
        for j=1%:size(xd,1)
%             Ddf(:,h,i,j)=ddf1(:,i,j).*Grad_y{j,i}(:,h+dim/2);
            Ddf(:,h,i,j)=ddf1(:,i,j).*Grad_y{j,i}(:,h);
        end
    end
end

grad_3=[];
for h=1:dim/2
    if type_num(h)==3
        uu(:,h)=1-cdf(dist_type{type_num(h)},-var_all(:,h)+2*X((h-1)*2+1),X((h-1)*2+1),X(h*2)); % 원래 seed 계산, Gumbel CDF
    else
        uu(:,h)=cdf(dist_type{type_num(h)},var_all(:,h),X((h-1)*2+1),X(h*2)); % 원래 seed 계산
    end
    grad_3=[grad_3 grad_chain3_ver02(uu(:,h),X((h-1)*2+1),X(h*2),type_num(h))];
end

% gg_1=[];grad_4=[];
% DDf=sum(Ddf,4);
% for i=1:1:size(v_exp_data,2)
% for h=1:dim/2
%     gg_1(i,[2*(h-1)+[1:2]])=[sum(repmat(DDf(:,h,i),1,2).*grad_3(:,(h-1)*2+1:h*2),1)];
% end
% end

gg_1=[];grad_4=[];
DDf=sum(Ddf,3);
% for i=1:1:size(v_exp_data,2)
for h=1:dim/2
    gg_1(:,[2*(h-1)+[1:2]])=[sum(repmat(DDf(:,h,1,1),1,2).*grad_3(:,(h-1)*2+1:h*2),1)];
end
% end

for h=1:dim/2
    grad_4=[grad_4 grad_chain4(X_((h-1)*2+1),X_(h*2),type_num(h))];
end

% gg=[];
% for i=1:1:size(v_exp_data,2)
% for h=1:dim/2
% gg(i,(h-1)*2+1:h*2)=[sum(gg_1(i,(h-1)*2+1:h*2).*grad_4(:,(h-1)*4+1:(h-1)*4+2),2) sum(gg_1(i,(h-1)*2+1:h*2).*grad_4(:,(h-1)*4+3:h*4),2)];
% end
% end

gg=[];
% for i=1:1:size(v_exp_data,2)
for h=1:dim/2
gg(:,(h-1)*2+1:h*2)=[sum(gg_1(:,(h-1)*2+1:h*2).*grad_4(:,(h-1)*4+1:(h-1)*4+2),2) sum(gg_1(:,(h-1)*2+1:h*2).*grad_4(:,(h-1)*4+3:h*4),2)];
end
% end

% gg_mult=sum(gg,1)/numel(pd_);
gg_mult=gg/numel(pd_);

if ff>=10^308
    gg=2*rand(dim,1)-1;
    gg_mult=gg.'/numel(pd_);
end
% disp('%%%%%%%%%%%%%%%%%%%% gradient evaluation %%%%%%%%%%%%%%%%%%%%')
time_analy=toc;
end
