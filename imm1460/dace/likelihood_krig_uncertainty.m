function obj=likelihood_krig_uncertainty(x)

for ii=1:1:size(x,1)
    X=x(ii,1:2);
disp(X)
load('exp_data004.mat')

[sim,beta]=model_krig(X);
sim=repmat(sim, [5,1]); %30

res(:,1)=(sim(:,1)-exp_data(:,1));
res(:,2)=(sim(:,2)-exp_data(:,2));
res(:,3)=(sim(:,3)-exp_data(:,3));

res=res.^2;

for i=1:1:3
res2(i)=(sum(res(:,i)))/beta(i);
% res2(i)=aa(i,ii).*(sum(res(:,i)))/beta(i);
end
% obj(ii,:)=exp(-5/2*log(beta(1)*beta(2)*beta(3)*beta(4)*beta(5)*beta(6)*beta(7)*beta(8))-0.5*(sum(res2)));
obj(ii,:)=(-5/2*log(2*pi)-5/2*log(beta(1)*beta(2)*beta(3))-0.5*(sum(res2)));
end

% obj=sum(res2);

