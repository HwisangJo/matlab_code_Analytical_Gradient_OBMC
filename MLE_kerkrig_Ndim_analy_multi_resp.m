%% Likelihood Function Metric
function [ff, gg]=MLE_kerkrig_Ndim_analy_multi_resp(x,func1,func2)

if nargout == 1 % gradient required
    [ff1]=func1(x);
    [ff2]=func2(x);
    ff=ff1+ff2;
elseif nargout>1
    [ff1, gg1]=func1(x);
    [ff2, gg2]=func2(x);
    ff=ff1+ff2;
    gg=gg1+gg2;
end

% if ff==inf
% disp('pause')
% end