function K = covBranchingProcess_3C(hyp, x, z, i)

%A three component branching process (where all branch from a latent process)

if nargin<2, K = '14'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists

xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

%Changepoint params.
location1  = hyp(1);
steepness1 = hyp(2);
location2  = hyp(3);
steepness2 = hyp(4);
location3  = hyp(5);
steepness3 = hyp(6);

covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}}; %Changepoint kernel
k1    = {'covMaterniso',3}; %Base kernel

if nargin<4 %Covariances  
    ind1 = find(x(:,2)==1);
    ind2 = find(x(:,2)==2);
    ind3 = find(x(:,2)==3);    
    if dg    
        K   = feval(k1{:},[hyp(13),hyp(14)],x(:,1));        
        K1   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x(ind1,1));
        K2   = feval(covN1{:}, [location2,steepness2,hyp(9),hyp(10)], x(ind2,1));
        K3   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind3,1));
        K(ind1,ind1) = K(ind1,ind1)+K1;
        K(ind2,ind2) = K(ind2,ind2)+K2;
        K(ind3,ind3) = K(ind3,ind3)+K3;                
        %K    = [K0(ind1,ind1)+K1,K0(ind1,ind2)  ,K0(ind1,ind3);
        %        K0(ind1,ind2)',K0(ind2,ind2)+K2 ,K0(ind2,ind3);
        %        K0(ind1,ind3)',K0(ind2,ind3)',K0(ind3,ind3)+K3];             
            
        K    = diag(K);
    else        
        if  xeqz   
        K    = feval(k1{:},[hyp(13),hyp(14)],x(:,1));
        covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}};
        %K1   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x(:,1));    
        %K2   = feval(covN1{:}, [location2,steepness2,hyp(9),hyp(10)], x(:,1));
        %K3   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(:,1));
        %K = [K0(ind1,ind1)+K1(ind1,ind1),K0(ind1,ind2),K0(ind1,ind3);
        %     K0(ind1,ind2)',K0(ind2,ind2)+K2(ind2,ind2),K0(ind2,ind3);
        %     K0(ind1,ind3)',K0(ind2,ind3)',K0(ind3,ind3)+K3(ind3,ind3)];
        K1   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x(ind1,1));    
        K2   = feval(covN1{:}, [location2,steepness2,hyp(9),hyp(10)], x(ind2,1));
        K3   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind3,1));
        %K = [K0(ind1,ind1)+K1,K0(ind1,ind2),K0(ind1,ind3);
        %     K0(ind1,ind2)',K0(ind2,ind2)+K2,K0(ind2,ind3);
        %     K0(ind1,ind3)',K0(ind2,ind3)',K0(ind3,ind3)+K3];         
        K(ind1,ind1) = K(ind1,ind1) + K1; 
        K(ind2,ind2) = K(ind2,ind2) + K2; 
        K(ind3,ind3) = K(ind3,ind3) + K3;         
        
        else %Cross covariances
        ind4 = find(z(:,2)==1);
        ind5 = find(z(:,2)==2);
        ind6 = find(z(:,2)==3);
        K    = feval(k1{:},[hyp(13),hyp(14)],x(:,1),z(:,1));
        %covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}};
        K1   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x(ind1,1),z(ind4,1));
        K2   = feval(covN1{:}, [location2,steepness2,hyp(9),hyp(10)], x(ind2,1),z(ind5,1));        
        K3   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind3,1),z(ind6,1));
        %K    = [K0(ind1,ind4)+K1,K0(ind1,ind5),K0(ind1,ind6);
        %        K0(ind2,ind4),K0(ind2,ind5)+K2,K0(ind2,ind6);
        %        K0(ind3,ind4),K0(ind3,ind5),K0(ind3,ind6)+K3];             
        try,K(ind1,ind4) = K(ind1,ind4)+K1;end
        try,K(ind2,ind5) = K(ind2,ind5)+K2;end     
        try,K(ind3,ind6) = K(ind3,ind6)+K3;end  
        %K1   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x(:,1),z(:,1));
        %K2   = feval(covN1{:}, [location2,steepness2,hyp(9),hyp(10)], x(:,1),z(:,1));        
        %K3   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(:,1),z(:,1));                        
        % K = [K0(ind1,ind4)+K1(ind1,ind4),K0(ind1,ind5),K0(ind1,ind6);
        %      K0(ind2,ind4),K0(ind2,ind5)+K2(ind2,ind5),K0(ind2,ind6);
        %      K0(ind3,ind4),K0(ind3,ind5),K0(ind3,ind6)+K3(ind3,ind6)];        
        end
end
    
else %Get derivatives (turn this into a loop for all labels)?
ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);
ind3 = find(x(:,2)==3);
  if i==1
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x(ind1,1),x(ind1,1),1);      
 elseif i==2
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x(ind1,1),x(ind1,1),2);          
  elseif i==3
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),1);      
   elseif i==4
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),2);  
 elseif i==5
    K = zeros(size(x,1),size(x,1));
    K(ind3,ind3)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind3,1),x(ind3,1),1);      
   elseif i==6
    K = zeros(size(x,1),size(x,1));
    K(ind3,ind3)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind3,1),x(ind3,1),2);        
  elseif i==7
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x(ind1,1),x(ind1,1),3);
  elseif i==8
    K = zeros(size(x,1),size(x,1));      
    K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x(ind1,1),x(ind1,1),4);    
  elseif i==9
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),3);
  elseif i==10
    K = zeros(size(x,1),size(x,1));      
    K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),4);    
  elseif i==11
    K = zeros(size(x,1),size(x,1));
    K(ind3,ind3)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind3,1),x(ind3,1),3);
  elseif i==12
    K = zeros(size(x,1),size(x,1));      
    K(ind3,ind3)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind3,1),x(ind3,1),4);    
  elseif i==13
     K = feval(k1{:},[hyp(13),hyp(14)],x(:,1),x(:,1),1);
  elseif i==14
     K = feval(k1{:},[hyp(13),hyp(14)],x(:,1),x(:,1),2);
  else
    error('Unknown hyperparameter')
  end
end