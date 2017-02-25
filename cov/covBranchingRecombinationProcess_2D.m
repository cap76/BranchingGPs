function K = covBranchingRecombinationProcess_2D(hyp, x, z, i)

% A two-component branching/recombination process covariance function

if nargin<2, K = '8'; return; end                 % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

locationb   = hyp(1); %Convergence time
%steepness1b = hyp(2); 
steepness2b = hyp(2); 
location   = hyp(3);  %Divergence time
%steepness1 = hyp(5);
steepness2 = hyp(4);

k1={'covMaterniso',3};

covN1 = {@covChangePointMultiD, {1, {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}}, @covZero}};

if nargin<4                                                        % covariances  
    ind1 = find(x(:,2)==1);
    ind2 = find(x(:,2)==2);
    
     if dg    
        K    = feval(k1{:},[hyp(5),hyp(6)],x(:,1));
        if isempty(ind2)==0        
        K3   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1));
        K(ind2,ind2) = K(ind2,ind2)+K3;   
        end
        K = diag(K);
     else        
        if  xeqz   
        K    = feval(k1{:},[hyp(5),hyp(6)],x(:,1));
         
        if isempty(ind2)==0
            K3   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1));            
            K(ind2,ind2) = K(ind2,ind2)+K3; 
        end
        
    else %Cross covariances
        ind3 = find(z(:,2)==1);
        ind4 = find(z(:,2)==2);        
        K    = feval(k1{:},[hyp(5),hyp(6)],x(:,1),z(:,1));
        
        if isempty(ind2)==0 & isempty(ind4)==0
            K3   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind4,1));            
            K(ind2,ind4) = K(ind2,ind4)+K3;
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
        K(ind2,ind4)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind4,1),1);                  
        else
        K = zeros(size(x,1),size(x,1));
        K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),1);          
        if dg,K = diag(K);end
        end
 
elseif i==2
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));K(ind2,ind4)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind4,1),2); 
        else
        K = zeros(size(x,1),size(x,1));
        K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),2);      
        if dg,K = diag(K);end
        end

    elseif i==3
         if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        K(ind2,ind4)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind4,1),3);                   
        else
        K = zeros(size(x,1),size(x,1));
        K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),3);                
        if dg,K = diag(K);end
        end

     elseif i==4
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));K(ind2,ind4)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind4,1),4); 
        else
        K = zeros(size(x,1),size(x,1));
        K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),4);              
        if dg,K = diag(K);end
        end
        
    elseif i==7
            if xeqz==0 & dg==0
            K = zeros(size(x,1),size(z,1)); K(ind2,ind4)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind4,1),5);
            else
            K = zeros(size(x,1),size(x,1));
            K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),5);
            if dg,K = diag(K);end
            end

    elseif i==8
            if xeqz==0 & dg==0
            K = zeros(size(x,1),size(z,1)); K(ind2,ind4)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind4,1),6);
            else
            K = zeros(size(x,1),size(x,1));
            K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),6);
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