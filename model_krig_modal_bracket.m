%% input variable 
function [Y,grad_y]=model_krig_modal_bracket(X,n_unknown)
load('dmodel_data_modal_bracket_200313.mat')

Y=zeros(size(X,1),1);
if nargout==1
    [Y]=predictor_mdfd_Ndim(X,dmodel1,n_unknown);
elseif nargout>1
    [Y grad_y]=predictor_mdfd_Ndim(X,dmodel1,n_unknown);
end
%%

%%
