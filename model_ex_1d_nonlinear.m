function [y,grad_y] = model_ex_1d_nonlinear(xx)

x1 = xx(:,1);

% y=exp(-1/(x.^2+4));
% y=log(x1.^2+4)+0.1*(x1-16).^3;
y=log(x1.^2+4)+0.1*cos(x1);

if nargout>1
%     grad_y=[exp(-1./(x.^2+4)).*(2*x)./(x.^2+4).^2];
%     grad_y=[2*x1./(x1.^2+4)+0.1*3*(x1-16).^2];
    grad_y=[2*x1./(x1.^2+4)-0.1*sin(x1)];
end
