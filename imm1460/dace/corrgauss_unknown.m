function  [r, dr] = corrgauss_unknown(theta, d)
%CORRGAUSS  Gaussian correlation function,
%
%           n
%   r_i = prod exp(-theta_j * d_ij^2) ,  i = 1,...,m
%          j=1
%
% If length(theta) = 1, then the model is isotropic:
% all  theta_j = theta .
%
% Call:    r = corrgauss(theta, d)
%          [r, dr] = corrgauss(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 

% hbn@imm.dtu.dk  
% Last update June 2, 2002

[m n] = size(d);  % number of differences and dimension of data
if  length(theta) == 1
  theta = repmat(theta,1,n);
elseif  length(theta) ~= n
  error(sprintf('Length of theta must be 1 or %d',n))
end

td = d.^2 .* repmat(-theta(:).',m,1);
r = exp(sum(td, 2));

unknown_n=1;
if  nargout > 1
%   dr = repmat(-2*theta(3:4),m,1) .* d(:,3:4) .* repmat(r,1,unknown_n);
    dr = repmat(-2*theta(n-unknown_n+1:n),m,1) .* d(:,n-unknown_n+1:n) .* repmat(r,1,unknown_n);
%   dr=[zeros(m,n-unknown_n) dr_];
end