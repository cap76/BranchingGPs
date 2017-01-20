function K = covBranchingProcess_2E(hyp, x, z, i)

%A two component branching process. Comprised of a base process and branch process.
%All Matern covariance functions.

%Should be fully compatible with gpml.

if nargin<2, K = '6'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

k1={'covMaterniso',3}; %Base kernel
 
location1  = hyp(1);
steepness1 = hyp(2);

covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}};

if nargin<4 %Covariances  
    ind1 = find(x(:,2)==1);
    ind2 = find(x(:,2)==2);
    
    if dg    
        K   = feval(k1{:},[hyp(5),hyp(6)],x(:,1));
        if isempty(ind2)==0
        K2   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1));
        K(ind2,ind2) = K(ind2,ind2)+K2;
        end
        K    = diag(K);
    else        
    if  xeqz   
    K    = feval(k1{:},[hyp(5),hyp(6)],x(:,1));
    
    if isempty(ind2)==0
    K2   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1));
    K(ind2,ind2) = K(ind2,ind2)+K2;
    end
    else %Cross covariances

        ind4 = find(z(:,2)==1);
        ind5 = find(z(:,2)==2);
         
        K    = feval(k1{:},[hyp(5),hyp(6)],x(:,1),z(:,1));
        if isempty(ind2)==0 & isempty(ind5)==0
        K2   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),z(ind5,1));
        K(ind2,ind5) = K(ind2,ind5)+K2;
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
    if isempty(ind2)==0 & isempty(ind4)==0,K(ind2,ind4)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),z(ind4,1),1);end        
    else
    K = zeros(size(x,1),size(x,1));
    if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),1);end      
    if dg,K = diag(K);end
    end
    
  elseif i==2
    
    if xeqz==0 & dg==0
    K = zeros(size(x,1),size(z,1));    
    if isempty(ind2)==0 & isempty(ind4)==0,K(ind2,ind4)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),z(ind4,1),2);end                   
    else      
    K = zeros(size(x,1),size(x,1));
    if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),2);end          
    if dg,K = diag(K);end
    end
    
    elseif i==3

    if xeqz==0 & dg==0
    K = zeros(size(x,1),size(z,1)); 
    if isempty(ind2)==0 & isempty(ind4)==0,K(ind2,ind4)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),z(ind4,1),3);end
    else     
    K = zeros(size(x,1),size(x,1));
    if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),3);end
    if dg,K = diag(K);end
    end
    
  elseif i==4
            
    if xeqz==0 & dg==0
    K = zeros(size(x,1),size(z,1)); 
    if isempty(ind2)==0 & isempty(ind4)==0,K(ind2,ind4)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),z(ind4,1),4);end           
    else
    K = zeros(size(x,1),size(x,1));      
    if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),4);end    
    if dg,K = diag(K);end
    end
  elseif i==5
     if xeqz==0 & dg==0 
     K = feval(k1{:},[hyp(5),hyp(6)],x(:,1),z(:,1),1);         
     else
     K = feval(k1{:},[hyp(5),hyp(6)],x(:,1),x(:,1),1);
     if dg,K = diag(K);end
     end
  elseif i==6
      
     if xeqz==0 & dg==0
     K = feval(k1{:},[hyp(5),hyp(6)],x(:,1),z(:,1),2);
     else
     K = feval(k1{:},[hyp(5),hyp(6)],x(:,1),x(:,1),2);
     if dg,K = diag(K);end
     end
     
  else
    error('Unknown hyperparameter')
  end
end