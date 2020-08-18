function  [f, df] = regpoly1_unknown(S)
%REGPOLY1  First order polynomial regression function
%
% Call:    f = regpoly1(S)
%          [f, df] = regpoly1(S)
%
% S : m*n matrix with design sites
% f = [1  s]
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update April 12, 2002

[m n] = size(S);
f = [ones(m,1)  S];

unknown_n=1;
if  nargout > 1
%   df = [zeros(1,2) eye(1)];
%     df = [zeros(n,1) eye(n)];
%     df(:,2:2+n-unknown_n-1)=zeros(n,n-unknown_n);
    df = zeros(unknown_n,n+1);
    df(end-unknown_n+1:end,end-unknown_n+1:end)=eye(unknown_n);
    df=repmat(df,m,1);
end