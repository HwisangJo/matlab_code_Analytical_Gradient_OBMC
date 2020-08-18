function [ab]=convert_param(mu,sd,num_t)
% num_t: distribution type
% 1: Normal, 2: Lognormal, 3: Gumbel, 5: Uniform, 6: Extreme

if num_t==1
    ab(1,:)=[mu sd];
elseif num_t==2
    ab(1,:)=[log(mu/sqrt(1+sd^2/mu^2)) sqrt(log(1+sd^2/mu^2))];
elseif num_t==3
    ab(1,:)=[mu-0.5772*sqrt(6)*sd/pi sqrt(6)*sd/pi];
elseif num_t==5
    ab(1,:)=[mu sd];
elseif num_t==6
    ab(1,:)=[mu+0.5772*sqrt(6)*sd/pi sqrt(6)*sd/pi];
end
