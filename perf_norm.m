%% Likelihood Function Metric
function [ff, gg]=perf_norm(x_n,obj,lbx,ubx)

x=(ubx-lbx).*x_n+lbx;

if nargout==1
ff=obj(x);
elseif nargout > 1 % gradient required
[ff gg_un]=obj(x);
gg=gg_un.*(ubx-lbx);
end
