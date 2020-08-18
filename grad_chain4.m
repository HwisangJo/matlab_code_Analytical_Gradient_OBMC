function g4=grad_chain4(mu,sd,num_t)
% num_t: distribution type
% 1: Normal, 2: Lognormal, 3: Gumbel, 5: Uniform, 6: Extreme

ab=convert_param(mu,sd,num_t);
a=ab(1);b=ab(2);

if num_t==1
    g4(1,:)=[1 0 0 1];
elseif num_t==2
    g4(1,:)=[1/mu*(1+sd^2/(mu^2+sd^2)) -sd^2/b/mu/(mu^2+sd^2) -sd/(mu^2+sd^2) sd/b/(mu^2+sd^2)];
elseif num_t==3
    g4(1,:)=[1 0 -0.5772*sqrt(6)/pi sqrt(6)/pi];
elseif num_t==5
    g4(1,:)=[1 0 0 1];
elseif num_t==6
    g4(1,:)=[1 0 0.5772*sqrt(6)/pi sqrt(6)/pi];
end
