function K = covConvergingProcess_2A(hyp, x, z, i)

%DUPLICATE
% A two component converging process with squared exponential covariance functions.

if nargin<2, K = '9'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode


ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);
 %keyboard

location = hyp(1);
steepness1 = hyp(2);
steepness2 = hyp(3);
%ell1 = exp(hyp(4));                                 % characteristic length scale
%sf21 = exp(2*hyp(5));                                           % signal variance
%ell2 = exp(hyp(4));                                 % characteristic length scale
%sf22 = exp(2*hyp(5));                                           % signal variance
%ell3 = exp(hyp(4));                                 % characteristic length scale
%sf23 = exp(2*hyp(5));                                           % signal variance

K1    = covSEiso([hyp(8),hyp(9)],x(:,1));
covN1 = {@covChangePointMultiD, {1, @covSEiso, @covZero}};
%covN2 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};

if nargin<4                                                        % covariances  

    %K = zeros(size(x,1),1);
    K2   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(:,1));
    K3   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(:,1));
    K    = [K2(ind1,ind1)+K1(ind1,ind1),K1(ind1,ind2);K1(ind2,ind1),K1(ind2,ind2)+K3(ind2,ind2)];           

    if dg
    K = diag(K);
    end
    
else                                                               % derivatives
  if i==1
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),1);      
    K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),x(ind2,1),1);          
  elseif i==2
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),2);          
  elseif i==3
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),x(ind2,1),2);  
  elseif i==4
    %K = zeros(size(x,1),size(x,1));
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),3);
  elseif i==5
    K = zeros(size(x,1),size(x,1));      
    %K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),4);
    K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),4);    
  elseif i==6
      %keyboard
      K = zeros(size(x,1),size(x,1));
     K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),x(ind2,1),3);
  elseif i==7
      K = zeros(size(x,1),size(x,1));
      K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),x(ind2,1),4);

  elseif i==8
     K = feval('covSEiso',[hyp(8),hyp(9)],x(:,1),x(:,1),1);
  elseif i==9
     K = feval('covSEiso',[hyp(8),hyp(9)],x(:,1),x(:,1),2);
  else
    error('Unknown hyperparameter')
  end
end