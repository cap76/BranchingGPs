function K = covBranchingProcess_2B(hyp, x, z, i)

% Branching process covariance function. A two-component branching processes 
% where both branches diverge from a latent process (SE-based).

%Should now fully function with gpml.

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
        
        if isempty(ind1)==0 
        K2   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1));
        K(ind1,ind1) = K(ind1,ind1) + K2;
        end
        
        if isempty(ind2)==0 
        K3   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1));
        K(ind2,ind2) = K(ind2,ind2) + K3;        
        end
        K = diag(K);
     else        
    if  xeqz   
    K    = covSEiso([hyp(8),hyp(9)],x(:,1));
    if isempty(ind1)==0
    K2   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1));
    K(ind1,ind1) = K(ind1,ind1) + K2;
    end
    if isempty(ind2)==0
    K3   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1));
    K(ind2,ind2) = K(ind2,ind2) + K3;
    end
    
    else %Cross covariances
        ind3 = find(z(:,2)==1);
        ind4 = find(z(:,2)==2);       
        K    = covSEiso([hyp(8),hyp(9)],x(:,1),z(:,1));

        if isempty(ind3)==0 & isempty(ind1)==0
        K2   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),z(ind3,1));
        K(ind1,ind3) = K(ind1,ind3) + K2;
        end
        
        if isempty(ind4)==0 & isempty(ind4)==0           
        K3   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),z(ind4,1));
        K(ind2,ind4) = K(ind2,ind4) + K3;
        end      
    end
    
end
    
else
    
ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);
if xeqz==0 & dg==0
    ind3 = find(z(:,2)==1);
    ind4 = find(z(:,2)==2);
end

% derivatives
  if i==1
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty(ind1)==0 & isempty(ind3)==0,K(ind1,ind3)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),z(ind3,1),1);end     
        if isempty(ind2)==0 & isempty(ind4)==0,K(ind2,ind4)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),z(ind4,1),1);end                  
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind1)==0,K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),1);end      
        if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),x(ind2,1),1);end          
        if dg,K = diag(K);end
        end

  elseif i==2
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        if isempty(ind1)==0 & isempty(ind3)==0,K(ind1,ind3)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),z(ind3,1),2);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind1)==0,K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),2);end          
        if dg,K = diag(K);end
        end

elseif i==3
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty(ind2)==0 & isempty(ind4)==0,K(ind2,ind4)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),z(ind4,1),2);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),x(ind2,1),2);end  
        if dg,K = diag(K);end
        end

    elseif i==4
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty(ind1)==0 & isempty(ind3)==0,K(ind1,ind3)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),z(ind3,1),3);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind1)==0,K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),3);end
        if dg,K = diag(K);end
        end

    elseif i==5
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty(ind1)==0 & isempty(ind3)==0,K(ind1,ind3)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),z(ind3,1),4);end 
        else
        K = zeros(size(x,1),size(x,1));      
        if isempty(ind1)==0,K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),4);end    
        if dg,K = diag(K);end
        end

    elseif i==6
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty(ind2)==0 & isempty(ind4)==0,K(ind2,ind4)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),z(ind4,1),3);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),x(ind2,1),3);end
        if dg,K = diag(K);end
        end

     elseif i==7
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));        
        if isempty(ind2)==0 & isempty(ind4)==0K(ind2,ind4)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),z(ind4,1),4);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(6),hyp(7)], x(ind2,1),x(ind2,1),4);end
        if dg,K = diag(K);end
        end

  elseif i==8
        if xeqz==0 & dg==0
        K = feval('covSEiso',[hyp(8),hyp(9)],x(:,1),z(:,1),1);
        else
        K = feval('covSEiso',[hyp(8),hyp(9)],x(:,1),x(:,1),1);
        if dg,K = diag(K);end
        end

     elseif i==9
        if xeqz==0 & dg==0
        K = feval('covSEiso',[hyp(8),hyp(9)],x(:,1),z(:,1),2);
        else
        K = feval('covSEiso',[hyp(8),hyp(9)],x(:,1),x(:,1),2);
        if dg,K = diag(K);end
        end

     else
    error('Unknown hyperparameter')
  end
end