function K = covBranchingProcess_2B(hyp, x, z, i)

% Branching process covariance function. Two-component branching processes 
% where both branches diverge from a latent process. 

if nargin<2, K = '9'; return; end                  % report number covBranchingProcess_2Bof parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

location = hyp(1);  
steepness1 = hyp(2);
steepness2 = hyp(3);
covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
  
if nargin<4                                                        % covariances  
    ind1 = find(x(:,2)==1);
    ind2 = find(x(:,2)==2);
    
     if dg    
        K    = covSEiso([hyp(8),hyp(9)],x(:,1));
        K2   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1));
        K3   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1));
        K(ind1,ind1) = K(ind1,ind1) + K2;
        K(ind2,ind2) = K(ind2,ind2) + K3;        
        K = diag(K);
     else        
    if  xeqz   
    K    = covSEiso([hyp(8),hyp(9)],x(:,1));
    K2   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1));
    K3   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1));
    K(ind1,ind1) = K(ind1,ind1) + K2;
    K(ind2,ind2) = K(ind2,ind2) + K3;            
    else %Cross covariances
        ind3 = find(z(:,2)==1);
        ind4 = find(z(:,2)==2);       
        K    = covSEiso([hyp(8),hyp(9)],x(:,1),z(:,1));

        if isempty(ind3)==0
        K2   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),z(ind3,1));
        K(ind1,ind3) = K(ind1,ind3) + K2;
        end
        
        if isempty(ind4)==0            
        K3   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),z(ind4,1));
        K(ind2,ind4) = K(ind2,ind4) + K3;
        end      
    end
    
end
    
else
    
ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);
  % derivatives
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