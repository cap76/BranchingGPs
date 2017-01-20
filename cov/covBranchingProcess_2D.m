function K = covBranchingProcess_2D(hyp, x, z, i)

%A two component branching process. Comprised of a latent (main) process from which
%two branches independetly branch.

if nargin<2, K = '10'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

k1={'covMaterniso',3}; %Base kernel
 
location1  = hyp(1);
steepness1 = hyp(2);
location2  = hyp(3);
steepness2 = hyp(4);

covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}};

if nargin<4 %Covariances  
    ind1 = find(x(:,2)==1);
    ind2 = find(x(:,2)==2);
    
    if dg    
        K   = feval(k1{:},[hyp(9),hyp(10)],x(:,1));
        K2   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1));
        K3   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1));
        K(ind1,ind1) = K(ind1,ind1)+K2;
        K(ind2,ind2) = K(ind2,ind2)+K3;
        K    = diag(K);
    else        
    if  xeqz   
    K    = feval(k1{:},[hyp(9),hyp(10)],x(:,1)); 
    
    if isempty(ind1)==0
    K2   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1));
    K(ind1,ind1) = K(ind1,ind1)+K2;
    end
    
    if isempty(ind2)==0
    K3   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1));           
    K(ind2,ind2) = K(ind2,ind2)+K3;
    end
    
    else %Cross covariances

    ind4 = find(z(:,2)==1);
    ind5 = find(z(:,2)==2);

         
    K    = feval(k1{:},[hyp(9),hyp(10)],x(:,1),z(:,1));
               
    if isempty(ind1)==0 & isempty(ind4)==0
    K2   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),z(ind4,1));        
    K(ind1,ind4) = K(ind1,ind4)+K2;
    end
            
    if isempty(ind2)==0 & isempty(ind5)==0
    K3   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind5,1));               
    K(ind2,ind5) = K(ind2,ind5)+K3;
    end
 
    end    
end
    
else
    
ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);

if xeqz==0 & dg==0
    ind4 = find(x(:,2)==1);
    ind5 = find(x(:,2)==2);
end

    %Note need to include the other options in here if we want to use FITC etc.
  % derivatives
  if i==1

        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));      K(ind1,ind4)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),z(ind4,1),1); 
        else
        K = zeros(size(x,1),size(x,1));
        K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),x(ind1,1),1);      
        if dg,K = diag(K);end
        end
        
        elseif i==2
            
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));      K(ind1,ind4)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),z(ind4,1),2);  
        else            
        K = zeros(size(x,1),size(x,1));
        K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),x(ind1,1),2);          
        if dg,K = diag(K);end
        end  
        
        elseif i==3
            
            
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));             K(ind2,ind5)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind5,1),1);       
        else            
        K = zeros(size(x,1),size(x,1));
        K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),1);      
        if dg,K = diag(K);end
        end  
        
        elseif i==4
            
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));      K(ind2,ind5)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind5,1),2);       
        else                        
        K = zeros(size(x,1),size(x,1));
        K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),2);          
        if dg,K = diag(K);end
        end
        
        elseif i==5

        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));  K(ind1,ind4)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),z(ind4,1),3);    
        else                        
        K = zeros(size(x,1),size(x,1));
        K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),x(ind1,1),3);
        if dg,K = diag(K);end
        end
          
        elseif i==6
           
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));  K(ind1,ind4)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),z(ind4,1),4);      
        else                        
        K = zeros(size(x,1),size(x,1));      
        K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),x(ind1,1),4);    
        if dg,K = diag(K);end
        end
        
        elseif i==7
        
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));    K(ind2,ind5)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind5,1),3);  
        else                   
        K = zeros(size(x,1),size(x,1));
        K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),3);
        if dg,K = diag(K);end
        end

        elseif i==8

        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));      K(ind2,ind5)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind5,1),4);
        else                   
        K = zeros(size(x,1),size(x,1));      
        K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),4);    
        if dg,K = diag(K);end
        end

        
        elseif i==9
            
        if xeqz==0 & dg==0
        K = feval(k1{:},[hyp(9),hyp(10)],x(:,1),z(:,1),1);
        else                                 
        K = feval(k1{:},[hyp(9),hyp(10)],x(:,1),x(:,1),1);
        if dg,K = diag(K);end
        end
     
  elseif i==10
        if xeqz==0 & dg==0
        K = feval(k1{:},[hyp(9),hyp(10)],x(:,1),z(:,1),2);
        else                               
        K = feval(k1{:},[hyp(9),hyp(10)],x(:,1),x(:,1),2);
        if dg,K = diag(K);end
        end
     
  else
    error('Unknown hyperparameter')
  end
end