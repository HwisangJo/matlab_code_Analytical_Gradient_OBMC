function [y,grad_y] = model_krig_mult_resp(x,y_num,func1,func2)
if y_num==1
    [y]=func1(x);
    if nargout>1
        [y,grad_y]=func1(x);
    end
elseif y_num==2
    [y]=func2(x);
    if nargout>1
        [y,grad_y]=func2(x);
    end
end
