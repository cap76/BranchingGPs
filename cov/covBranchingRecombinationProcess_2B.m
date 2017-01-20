function K = covBranchingRecombinationProcess_2B(hyp, x, z, i)

%THIS IS JUNK ATM!
% A 2-component branching/recombination process covariance function
%This will be a two-component branch/recomb covariance function with the branch
%recomb independent.

if nargin<2, K = '12'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

locationb   = hyp(1); %Convergence time
steepness1b = hyp(2); 
steepness2b = hyp(3); 
location    = hyp(4);
steepness1  = hyp(5);
steepness2  = hyp(6);

covN1 = {@covChangePointMultiD, {1, {@covChangePointMultiD, {1, @covZero, @covSEiso}}, @covZero}};

ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);  

if nargin<4                                                        % covariances  
    
     if dg    
        K1   = covSEiso([hyp(11),hyp(12)],x(:,1));
        K2   = feval(covN1{:}, [locationb,steepness1b,location,steepness1,hyp(7),hyp(8)], x(:,1));
        K3   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(9),hyp(10)], x(:,1));
        K    = [K1(ind1,ind1),K1(ind1,ind2);K1(ind2,ind1),K1(ind2,ind2)+K3(ind2,ind2)];
        K    = diag(K);
     else        
        if  xeqz   
        K1    = covSEiso([hyp(11),hyp(12)],x(:,1));
        K2   = feval(covN1{:}, [locationb,steepness1b,location,steepness1,hyp(7),hyp(8)], x(:,1));
        K3   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(9),hyp(10)], x(:,1));
        K    = [K1(ind1,ind1),K1(ind1,ind2);K1(ind2,ind1),K1(ind2,ind2)+K3(ind2,ind2)];           
        else %Cross covariances
        ind3 = find(z(:,2)==1);
        ind4 = find(z(:,2)==2);        
        K1    = covSEiso([hyp(11),hyp(12)],x(:,1),z(:,1));
        K2   = feval(covN1{:}, [locationb,steepness1b,location,steepness1,hyp(7),hyp(8)], x(:,1),z(:,1));
        K3   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(9),hyp(10)], x(:,1),z(:,1));        
        K = [K1(ind1,ind3),K1(ind1,ind4); K1(ind2,ind3),K1(ind2,ind4)+K3(ind2,ind4)];        
        end    
end
    
else
    
if xeqz==0 & dg==0
    ind3 = find(z(:,2)==1);
    ind4 = find(z(:,2)==2);
end

  % derivatives
  if i==1
      
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),1);          

  elseif i==2
    K = zeros(size(x,1),size(x,1));
    %K(ind1,ind1)   = feval(covN1{:}, [locationb,steepness1b,location,steepness1,hyp(7),hyp(8)], x(ind1,1),x(ind1,1),2);          
  elseif i==3
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),2);      
  elseif i==4
    K = zeros(size(x,1),size(x,1));
    %K(ind1,ind1)   = feval(covN1{:}, [locationb,steepness1b,location,steepness1,hyp(7),hyp(8)], x(ind1,1),x(ind1,1),3);      
    K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),3);                
  elseif i==5
     K = zeros(size(x,1),size(x,1));
     %K(ind1,ind1)   = feval(covN1{:}, [locationb,steepness1b,location,steepness1,hyp(7),hyp(8)], x(ind1,1),x(ind1,1),4);          
  elseif i==6
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),4);          
    
  elseif i==7
    %K = zeros(size(x,1),size(x,1));
    K = zeros(size(x,1),size(x,1));
    %K(ind1,ind1)   = feval(covN1{:}, [locationb,steepness1b,location,steepness1,hyp(7),hyp(8)], x(ind1,1),x(ind1,1),5);
  elseif i==8
    K = zeros(size(x,1),size(x,1));      
    %K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),4);
    %K(ind1,ind1)   = feval(covN1{:}, [locationb,steepness2b,location,steepness1,hyp(7),hyp(8)], x(ind1,1),x(ind1,1),6);    
  elseif i==9
      %keyboard
      K = zeros(size(x,1),size(x,1));
     K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness1b,location,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),5);
  elseif i==10
      K = zeros(size(x,1),size(x,1));
      K(ind2,ind2)   = feval(covN1{:}, [locationb,steepness2b,location,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),6);
  elseif i==11
     K = feval('covSEiso',[hyp(11),hyp(12)],x(:,1),x(:,1),1);
  elseif i==12
     K = feval('covSEiso',[hyp(11),hyp(12)],x(:,1),x(:,1),2);
  else
    error('Unknown hyperparameter')
  end
end