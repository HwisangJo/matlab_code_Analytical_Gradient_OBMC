%% input variable 
function [Y,grad_y_]=model_krig_modal_bracket2(X,n_unknown)
load('dmodel_data_modal_bracket_200313.mat')

Y=zeros(size(X,1),1);
if nargout==1
    [Y]=predictor_mdfd_Ndim(X(:,1),dmodel2,n_unknown);
elseif nargout>1
    [Y grad_y]=predictor_mdfd_Ndim(X(:,1),dmodel2,n_unknown);
    grad_y_=[grad_y zeros(size(grad_y,1),size(grad_y,2))];
end



%%

%%
