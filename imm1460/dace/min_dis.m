%% minimum distance
function y=min_dis(XX)
load('seq_samp_prep.mat');
x=x_B;
dis=[];
for i=1:length(x)
dis(i)=sum((x(i,:)-XX).^2);
end
y=min(dis);
%%
