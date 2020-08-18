%% input variable 
function [y,beta]=model_krig(X)
load('krig_model.mat')

[result1,dy1,mse1]=predictor(X,dmodel1);
[result2,dy2,mse2]=predictor(X,dmodel2);
[result3,dy3,mse3]=predictor(X,dmodel3);

y=[result1 result2 result3];
beta=[mse1,mse2,mse3];
dyy=[dy1,dy2,dy3];
%%
