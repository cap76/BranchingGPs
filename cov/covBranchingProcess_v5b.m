function K = covBranchingProcess_v5b(hyp, x, z, i)

%A two component branching process.

if nargin<2, K = '6'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

k1={'covMaterniso',3}; %Base kernel
 
location1  = hyp(1);
steepness1 = hyp(2);
covN1 = {@covChangePointMultiD, {1, {'covMaterniso',3},@covZero}};
%covN2 = {@covSum,{covSEiso,{@covChangePointMultiD, {1, @covZero, @covSEiso}}}};
%covN3 = {@covSum,{covSEiso,{@covSum,{{@covChangePointMultiD, {1, @covZero, @covSEiso}},{@covChangePointMultiD, {1, @covZero, @covSEiso}}}}}};

if nargin<4 %Covariances  
    ind1 = find(x(:,2)==1);
    ind2 = find(x(:,2)==2);
    
    if dg    
        K1   = feval(k1{:},[hyp(5),hyp(6)],x(:,1));
        K2   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(:,1));
        K    = [K1(ind1,ind1),K1(ind1,ind2);
                K1(ind1,ind2)',K1(ind2,ind2)+K2(ind2,ind2)];        
        K    = diag(K);
    else        
    if  xeqz   
    K1    = feval(k1{:},[hyp(5),hyp(6)],x(:,1));
    covN1 = {@covChangePointMultiD, {1, {'covMaterniso',3}, @covZero}};
    K2   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(:,1));
    K = [K1(ind1,ind1),K1(ind1,ind2);
         K1(ind1,ind2)',K1(ind2,ind2)+K2(ind2,ind2)];
    else %Cross covariances
        ind4 = find(z(:,2)==1);
        ind5 = find(z(:,2)==2);

        
        K1    = feval(k1{:},[hyp(5),hyp(6)],x(:,1),z(:,1));
        covN1 = {@covChangePointMultiD, {1, {'covMaterniso',3}, @covZero}};
        %covN2 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
        %K = zeros(size(x,1),1);
        K2   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(:,1),z(:,1));
        
         K = [K1(ind1,ind4),K1(ind1,ind5);            
             K1(ind2,ind4),K1(ind2,ind5)+K2(ind2,ind5)];
        
    end
    
end
    
else
    
ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);
  % derivatives
  if i==1
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),1);      
    %K(ind3,ind3)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind3,1),x(ind3,1),1);      
    %K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),1);          
 elseif i==2
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),2);          
    %K(ind3,ind3)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind3,1),x(ind3,1),2);              
  elseif i==3
    %K = zeros(size(x,1),size(x,1));
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),3);
    %K(ind3,ind3)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind3,1),x(ind3,1),3);    
  elseif i==4
    K = zeros(size(x,1),size(x,1));      
    %K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),4);
    K(ind2,ind2)   = feval(covN1{:}, [location1,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),4);    
    %K(ind3,ind3)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind3,1),x(ind3,1),4);        
  elseif i==5
     K = feval(k1{:},[hyp(5),hyp(6)],x(:,1),x(:,1),1);
  elseif i==6
     K = feval(k1{:},[hyp(5),hyp(6)],x(:,1),x(:,1),2);
  else
    error('Unknown hyperparameter')
  end
end