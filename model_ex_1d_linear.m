function [y,grad_y] = model_ex_1d_linear(xx)

x1 = xx(:,1);

y=1.25*x1;

if nargout>1
    grad_y=1.25;
end
