%% CVE Kriging model
function y=cve_krig(XX)

load('seq_samp_prep3.mat');

result1=predictor(XX,dmodel1);
result2=predictor(XX,dmodel2);
result3=predictor(XX,dmodel3);

y=[result1];