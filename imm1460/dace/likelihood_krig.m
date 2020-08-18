function obj=likelihood_krig(x)

for ii=1:1:size(x,1)
    X=x(ii,1:2);
disp(X)
load('exp_data004.mat')

% exp_data=exp_data004;

% sensitiv=[0	 0	6.94553E-13	8.84763E-14	0.000876576	6.45978E-15	9.73969E-13	2.51449E-14
%           0.998207455  	0.000908806	1.01207E-16	1.17022E-13	7.16315E-06	6.30296E-14	1.10469E-14	4.8279E-13];
% sensitiv=1*ones(8,2);
% aa(:,ii)=sensitiv*delta';
%  kk(ii)=sum(aa(:,ii));
%  aa(:,ii)=aa(:,ii)./kk(ii);
% beta(1)=var(exp_data(:,1));
% beta(2)=var(exp_data(:,2));
% beta(3)=var(exp_data(:,3));
% beta(4)=var(exp_data(:,4));
% beta(5)=var(exp_data(:,5));
% beta(6)=var(exp_data(:,6));
% beta(7)=var(exp_data(:,7));
% beta(8)=var(exp_data(:,8));

[sim, bet]=model_krig(X);
sim=repmat(sim, [5,1]); %30

res(:,1)=(sim(:,1)-exp_data(:,1));
res(:,2)=(sim(:,2)-exp_data(:,2));
res(:,3)=(sim(:,3)-exp_data(:,3));

% 
% beta(1)=var(res(:,1))*5/3;
% beta(2)=var(res(:,2))*5/3;
% beta(3)=var(res(:,3))*5/2;
% beta(4)=1;%var(res(:,4))*5/2;
% beta(5)=1;%var(res(:,5))*5/2;
% beta(6)=1;%var(res(:,6))*5/2;
% beta(7)=1;%var(res(:,7))*5/2;
% beta(8)=1;%var(res(:,8))*5/2;
% beta=[7.916666667   1.319444444   3.15206E-06   0.000169   0.038173605   1.91406E-05   3.55313E-06   0.000169];


% for i=1:1:8;
%     beta(i)=0.1;
% end

% beta(1)=sum(res(:,1).^2)/3;
% beta(2)=sum(res(:,2).^2)/3.^0.5;
% beta(3)=sum(res(:,3).^2)/3.^0.5;
% beta(4)=sum(res(:,4).^2)/3.^0.5;
% beta(5)=sum(res(:,5).^2)/3.^0.5;
% beta(6)=sum(res(:,6).^2)/3.^0.5;
% beta(7)=sum(res(:,7).^2)/3.^0.5;
% beta(8)=sum(res(:,8).^2)/3.^0.5;


% 
beta(1)=(sum((res(:,1)-mean(res(:,1))).^2))/3;
beta(2)=(sum((res(:,2)-mean(res(:,2))).^2))/3;
beta(3)=(sum((res(:,3)-mean(res(:,3))).^2))/3;


% beta=[1.666675360563723   1.666928542951552   1.785044035751308   1.695503290358453   1.668218791907860   1.634166140430496   1.602494097656588   1.643144004964273];


% % a=3/0.072;
% % beta(1)=a*8.842313472806235e-06;
% % beta(2)=a*4.925755060583665e-05;
% % beta(3)=a*8.459149574630842e-06;
% % beta(4)=a*1.489846252581402e-05;
% % beta(5)=a*6.198771800688408e-07; %1;%
% % beta(6)=a*4.444426153031078e-05;
% % beta(7)=a*1.872256026319736e-06;
% % beta(8)=a*6.615288375928616e-06; %1;%

res=res.^2;

for i=1:1:3
res2(i)=(sum(res(:,i)))/beta(i);
% res2(i)=aa(i,ii).*(sum(res(:,i)))/beta(i);
end
% obj(ii,:)=exp(-5/2*log(beta(1)*beta(2)*beta(3)*beta(4)*beta(5)*beta(6)*beta(7)*beta(8))-0.5*(sum(res2)));
obj(ii,:)=(-5/2*log(2*pi)-5/2*log(beta(1)*beta(2)*beta(3))-0.5*(sum(res2)));
end

% obj=sum(res2);

