% function  [f, df] = regpoly2_unknown(S)
function  [f, ddf] = regpoly2_unknown(S)
%REGPOLY2  Second order polynomial regression function
% Call:    f = regpoly2(S)
%          [f, df] = regpoly2(S)
%
% S : m*n matrix with design sites
% f =  [1 S S(:,1)*S S(:,2)S(:,2:n) ... S(:,n)^2]
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update September 4, 2002

[m n] = size(S);
nn = (n+1)*(n+2)/2;  % Number of columns in f  
% Compute  f
f = [ones(m,1) S zeros(m,nn-n-1)];
j = n+1;   q = n;
for  k = 1 : n
  f(:,j+(1:q)) = repmat(S(:,k),1,q) .* S(:,k:n);
  j = j+q;   q = q-1;
end

unknown_n=1;
if  nargout > 1
  df = [zeros(n,1)  eye(n)  zeros(n,nn-n-1)];
  df(:,2:2+n-unknown_n-1)=zeros(n,n-unknown_n);
  df=repmat(df,m,1);
  j = n+1;   q = n; 
  for  k = 1 : n
%       if k>=n-unknown_n+1 && k<=n
    df(k:n:n*(m-1)+k,j+(1:q)) = [2*S(:,k) S(:,k+1:n)];
%     end
    for i = 1 : n-k
%         if k>=n-unknown_n+1 && k<=n
        df([k:n:n*(m-1)+k]+i,j+1+i) = S(:,k);
%         end
    end
    j = j+q;   q = q-1;
  end
  
  ddf=[];
  for k=n-unknown_n+1:n
  ddf=[ddf;df([k:n:n*(m-1)+k],:)];
  end
end 