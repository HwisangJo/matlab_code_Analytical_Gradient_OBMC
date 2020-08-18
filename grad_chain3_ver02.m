function g3=grad_chain3_ver02(u,a,b,num_t)
% num_t: distribution type
% 1: Normal, 2: Lognormal, 3: Gumbel, 4: Weibull, 5: Uniform, 6: Extreme

if num_t==1
    g3=[ones(size(u,1),1) norminv(u)];
elseif num_t==2
    g3=[exp(a+b*norminv(u)) exp(a+b*norminv(u)).*norminv(u)];
elseif num_t==3
    g3=[ones(size(u,1),1) -log(-log(u))];
elseif num_t==4
    g3=[(-log(1-u)).^(1/b) -a/b^2*(-log(1-u)).^(1/b).*log(-log(1-u))];
elseif num_t==5
    g3=[ones(size(u,1),1) sqrt(3)*(2*u-1)];
elseif num_t==6
    g3=[ones(size(u,1),1) log(-log(1-u))];
end