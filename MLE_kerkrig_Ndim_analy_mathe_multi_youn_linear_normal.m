%% Likelihood Function Metric
function [ff, gg_mult,time_analy]=MLE_kerkrig_Ndim_analy_mathe_multi_youn_linear_normal(X_,file_exp_name,perf_func)
% ex) type_num=[1 3 5 2]; % for 4 unknown variables
% ex) X_=[mu1 sigma1 ... muN sigmaM]; % for N unknown variables
% 1: Normal, 2: Lognormal, 3: Gumbel, 5: Uniform, 6: Extreme
% perf_func=@(x,n_unknown) model_krig_10bar_ver02(x,n_unknown);
load(file_exp_name)

dim=size(X_,2)/2;
mean_x=X_(2*[1:dim]-1);
% std_x=X_(2*[1:dim]);
var_x=X_(2*[1:dim]);

% Uncertainty propagation(KDE) & Likelihood Function Metric
for i=1:1:size(v_exp_data,2) %y
%     Xd=xd(:,i);
    Xd=[];
%     parfor j=1%:size(xd,1)
    for j=1%:size(xd,1)
%         [y,grad_y]=perf_func([repmat(Xd(j,:),n,1), var_all],i);
        if size(v_exp_data,2)>1
        [y,grad_y]=perf_func(mean_x,i);
        elseif size(v_exp_data,2)==1
        [y,grad_y]=perf_func(mean_x);
        grad_y=grad_y(1);
        end
        simm(:,j,i)=y;
        Grad_y{j,i}=grad_y;

%         mean_y=mean_x*grad_y'+y;
        mean_y=y;
        var_y=var_x*(grad_y.^2)';
    end
end

N=size(v_exp_data,3);
ff=sum((v_exp_data-mean_y).^2)/(2*var_y)+log(var_y)*N/2+N*log(2*pi)/2;

if ff>10^308
    ff=10^308;
end
% disp('%%%%%%%%%%%%%%%%%%%% function evaluation %%%%%%%%%%%%%%%%%%%%')
if nargout > 1 % gradient required
%     tic
gg_mean=-grad_y*sum(v_exp_data-mean_y)/var_y;
gg_std=-grad_y.^2*sum((v_exp_data-mean_y).^2)/(2*var_y^2)+N*grad_y.^2/(2*var_y);
% gg_std=abs(grad_y)*(-sum((v_exp_data-mean_y).^2)/var_y^(3/2)+N/sqrt(var_y));

gg_mult=zeros(1,dim*2);
gg_mult(2*[1:dim]-1)=gg_mean;
gg_mult(2*[1:dim])=gg_std;

if ff>=10^308
    gg=2*rand(dim,1)-1;
    gg_mult=gg/numel(pd_);
end
% disp('%%%%%%%%%%%%%%%%%%%% gradient evaluation %%%%%%%%%%%%%%%%%%%%')
% time_analy=toc;
end
