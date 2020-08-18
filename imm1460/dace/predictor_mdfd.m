function  [y, or1, or2, dmse] = predictor_mdfd(x, dmodel)
%PREDICTOR  Predictor for y(x) using the given DACE model.
%
% Call:   y = predictor(x, dmodel)
%         [y, or] = predictor(x, dmodel)
%         [y, dy, mse] = predictor(x, dmodel) 
%         [y, dy, mse, dmse] = predictor(x, dmodel) 
%
% Input
% x      : trial design sites with n dimensions.  
%          For mx trial sites x:
%          If mx = 1, then both a row and a column vector is accepted,
%          otherwise, x must be an mx*n matrix with the sites stored
%          rowwise.
% dmodel : Struct with DACE model; see DACEFIT
%
% Output
% y    : predicted response at x.
% or   : If mx = 1, then or = gradient vector/Jacobian matrix of predictor
%        otherwise, or is an vector with mx rows containing the estimated
%                   mean squared error of the predictor
% Three or four results are allowed only when mx = 1,
% dy   : Gradient of predictor; column vector with  n elements
% mse  : Estimated mean squared error of the predictor;
% dmse : Gradient vector/Jacobian matrix of mse

% hbn@imm.dtu.dk
% Last update August 26, 2002
 
  or1 = NaN;   or2 = NaN;  dmse = NaN;  % Default return values
  if  isnan(dmodel.beta)
    y = NaN;   
    error('DMODEL has not been found')
  end

  [m n] = size(dmodel.S);  % number of design sites and number of dimensions
  sx = size(x);            % number of trial sites and their dimension
  if  min(sx) == 1 & n > 1 % Single trial point 
    nx = max(sx);
    if  nx == n 
      mx = 1;  x = x(:).';
    end
  else
    mx = sx(1);  nx = sx(2);
  end
  if  nx ~= n
    error(sprintf('Dimension of trial sites should be %d',n))
  end
  
  % Normalize trial sites  
  x = (x - repmat(dmodel.Ssc(1,:),mx,1)) ./ repmat(dmodel.Ssc(2,:),mx,1);
  q = size(dmodel.Ysc,2);  % number of response functions
  y = zeros(mx,q);         % initialize result
  
  if  mx == 1  % one site only
    dx = repmat(x,m,1) - dmodel.S;  % distances to design sites
    if  nargout > 1                 % gradient/Jacobian wanted
      [f df] = feval(dmodel.regr, x);
      [r dr] = feval(dmodel.corr, dmodel.theta, dx);
      % Scaled Jacobian
      drr=[];
      for ik=1:mx
        drr=[drr dr((ik-1)*m+1:ik*m,:)];
      end
      
      dy = (repmat(df,mx,1) * dmodel.beta).' + dmodel.gamma * drr;
      
      dyy=[];
      for ik=1:mx
        dyy=[dyy; dy((ik-1)*n+1:ik*n)];
      end
      
      % Unscaled Jacobian
      or1 = dyy .* repmat(dmodel.Ysc(2, :)', 1, nx) ./ repmat(dmodel.Ssc(2,:), q, 1);
      if q == 1
        % Gradient as a column vector
        or1 = or1';
      end
      if  nargout > 2  % MSE wanted
        
        rt = dmodel.C \ r;
        u = dmodel.Ft.' * rt - f.';
        v = dmodel.G \ u;
        or2 = repmat(dmodel.sigma2,mx,1) .* repmat((1 + sum(v.^2) - sum(rt.^2))',1,q);
        
        if  nargout > 3  % gradient/Jacobian of MSE wanted
          % Scaled gradient as a row vector
          Gv = dmodel.G' \ v;
          g = (dmodel.Ft * Gv - rt)' * (dmodel.C \ dr) - (df * Gv)';
          % Unscaled Jacobian
          dmse = repmat(2 * dmodel.sigma2',1,nx) .* repmat(g ./ dmodel.Ssc(2,:),q,1);
          if q == 1
            % Gradient as a column vector
            dmse = dmse';
          end
        end
        
      end
      
    else  % predictor only
      f = feval(dmodel.regr, x);
      r = feval(dmodel.corr, dmodel.theta, dx);
    end
    
    % Scaled predictor
    sy = f * dmodel.beta + (dmodel.gamma*r).';
    % Predictor
    y = (dmodel.Ysc(1,:) + dmodel.Ysc(2,:) .* sy)';
    
  else  % several trial sites
    % Get distances to design sites  

    % Get regression function and correlation
    if nargout==1
        dx = zeros(mx*m,n);  kk = 1:m;
        for  k = 1 : mx
            dx(kk,:) = repmat(x(k,:),m,1) - dmodel.S;
            kk = kk + m;
        end
        f = feval(dmodel.regr, x);
        r = feval(dmodel.corr, dmodel.theta, dx);
        r = reshape(r, m, mx);
    elseif nargout>1
%         dx = zeros(mx*m,n);  kk = 1:mx;
%         for  k = 1 : m
%             dx(kk,:) = x - dmodel.S(k,:);
%             kk = kk + mx;
%         end
%         [f, df] = feval(dmodel.regr, x);
%         [r, dr] = feval(dmodel.corr, dmodel.theta, dx);
%         r = reshape(r, mx, m)';
        
        dx = zeros(mx*m,n);  kk = 1:m;
        for  k = 1 : mx
            dx(kk,:) = repmat(x(k,:),m,1) - dmodel.S;
            kk = kk + m;
        end
%         [f, df] = feval(dmodel.regr, x);
%         [r, dr] = feval(dmodel.corr, dmodel.theta, dx);

        [f, df] = regpoly1_unknown(x);
%         [f, df] = regpoly2_unknown(x);
        [r, dr] = corrgauss_unknown(dmodel.theta, dx);
        r = reshape(r, m, mx);
    end
    
    % Scaled predictor 
    sy = f * dmodel.beta + (dmodel.gamma * r).';
    % Predictor
    y = repmat(dmodel.Ysc(1,:),mx,1) + repmat(dmodel.Ysc(2,:),mx,1) .* sy;
    
    if  nargout > 1   % MSE wanted
%       rt = dmodel.C \ r;
%       u = dmodel.G \ (dmodel.Ft.' * rt - f.');
%       or2 = repmat(dmodel.sigma2,mx,1) .* repmat((1 + colsum(u.^2) - colsum(rt.^2))',1,q);

%       [f df] = feval(dmodel.regr, x);
%       [r, dr] = feval(dmodel.corr, dmodel.theta, dx);
%       r = reshape(r, m, mx);
%       dr = corrgauss_dr(dmodel.theta, dx);
      
      % Scaled Jacobian
%       drr=[];
%       for ik=1:mx
%         drr=[drr dr((ik-1)*m+1:ik*m,:)];
%       end

%       drr=zeros(m,mx*n);
      drr=zeros(m,mx*2);
%       tic
%       for ik=1:2 %n-unknown_n:n
    for ik=1%:2 %n-unknown_n:n
%         drr(:,[ik:2:(mx-1)*2+ik])=reshape(dr(:,ik),m,mx);
        drr=reshape(dr(:,ik),m,mx);
      end
      
%       for ik=1:n %n-unknown_n:n
%         drr(:,[ik:n:(mx-1)*n+ik])=reshape(dr(:,ik),m,mx);
%       end
      
%       toc

% tic
%        drr=reshape(dr',n*mx,m)';
% toc
       
%       dr_=reshape(dr,m,mx*n);drr=[];
%       for ik=1:mx
%         drr=[drr dr_(:,[ik:mx:mx*(n-1)+ik])];
%       end

      dy = (repmat(df,mx,1) * dmodel.beta).' + dmodel.gamma * drr;
%       dy = (df * dmodel.beta).' + dmodel.gamma * drr;
      
%       dyy=[];
% tic
%       for ik=1:n
%         dyy=[dyy dy(ik:n:(mx-1)*n+ik)'];
%       end
% toc

% tic
%       dyy=reshape(dy',n,mx)';
        dyy=reshape(dy',1,mx)';
% toc
      
      % Unscaled Jacobian
      or1 = dyy .* repmat(dmodel.Ysc(2, :)', 1, nx-1) ./ repmat(dmodel.Ssc(2,2), q, 1);
      if q == 1
        % Gradient as a column vector
%         or1 = or1';
      end
      if  nargout > 2
        disp('WARNING from PREDICTOR.  Only  y  and  or1=mse  are computed')
      end
    end
    
  end % of several sites
  
% >>>>>>>>>>>>>>>>   Auxiliary function  ====================

function  s = colsum(x)
% Columnwise sum of elements in  x
if  size(x,1) == 1,  s = x; 
else,                s = sum(x);  end