function K = covBranchingProcess_3A(hyp, x, z, i)

%A three component branching process with a main branch with two processess 
%branching in a reciprocal manner.

if nargin<2, K = '10'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

location1  = hyp(1);
steepness1 = hyp(2);
location2  = hyp(3);
steepness2 = hyp(4);
covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
%covN2 = {@covSum,{covSEiso,{@covChangePointMultiD, {1, @covZero, @covSEiso}}}};
%covN3 = {@covSum,{covSEiso,{@covSum,{{@covChangePointMultiD, {1, @covZero, @covSEiso}},{@covChangePointMultiD, {1, @covZero, @covSEiso}}}}}};

if nargin<4 %Covariances  
    ind1 = find(x(:,2)==1);
    ind2 = find(x(:,2)==2);
    ind3 = find(x(:,2)==3);
    
    if dg    
        K   = covSEiso([hyp(9),hyp(10)],x(:,1));
        K2   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x([ind2;ind3],1));
        K3   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind3,1));
        K([ind2;ind3],[ind2;ind3]) = K([ind2;ind3],[ind2;ind3])+K2;
        K(ind3,ind3) = K(ind3,ind3)+K3;        
        
        %K    = [K1(ind1,ind1),K1(ind1,ind2)               ,   K1(ind1,ind3);
        %        K1(ind1,ind2)',K1(ind2,ind2)+K2(ind2,ind2) , K1(ind2,ind3)+K2(ind2,ind3);
        %        K1(ind1,ind3)',K1(ind2,ind3)'+K2(ind2,ind3)' ,K1(ind3,ind3)+K2(ind3,ind3)+K3(ind3,ind3)];        
        K    = diag(K);
    else        
    if  xeqz   
    K    = covSEiso([hyp(9),hyp(10)],x(:,1));
    covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
    K2   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x([ind2;ind3],1));
    K3   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind3,1));
    K([ind2;ind3],[ind2;ind3]) = K([ind2;ind3],[ind2;ind3])+K2;
    K(ind3,ind3) = K(ind3,ind3)+K3;    
    %K = [K1(ind1,ind1),K1(ind1,ind2),K1(ind1,ind3);
    %         K1(ind1,ind2)',K1(ind2,ind2)+K2(ind2,ind2),K1(ind2,ind3)+K2(ind2,ind3);
    %         K1(ind1,ind3)',K1(ind2,ind3)'+K2(ind2,ind3)',K1(ind3,ind3)+K2(ind3,ind3)+K3(ind3,ind3)];
    else %Cross covariances
        ind4 = find(z(:,2)==1);
        ind5 = find(z(:,2)==2);
        ind6 = find(z(:,2)==3);

        
        K    = covSEiso([hyp(9),hyp(10)],x(:,1),z(:,1));
        covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
        %covN2 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
        %K = zeros(size(x,1),1);
        K2   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x([ind2;ind3],1),z([ind5;ind6],1));
        K3   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind3,1),z(ind6,1));
        K([ind2;ind3],[ind5;ind6]) = K([ind2;ind3],[ind5;ind6])+K2;
        K(ind3,ind6) = K(ind3,ind6)+K3;        
        % K = [K1(ind1,ind4),K1(ind1,ind5),K1(ind1,ind6);
        %     K1(ind2,ind4),K1(ind2,ind5)+K2(ind2,ind5),K1(ind2,ind6)+K2(ind2,ind6);
        %     K1(ind3,ind4),K1(ind3,ind5)+K2(ind3,ind5),K1(ind3,ind6)+K2(ind3,ind6)+K3(ind3,ind6)];
        
    end
    
end
    
else
    
ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);
ind3 = find(x(:,2)==3);
  % derivatives
  if i==1
    K = zeros(size(x,1),size(x,1));
    K([ind2;ind3],[ind2;ind3])   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x([ind2;ind3],1),x([ind2;ind3],1),1);      
    %K(ind3,ind3)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind3,1),x(ind3,1),1);      
    %K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),1);          
 elseif i==2
    K = zeros(size(x,1),size(x,1));
    K([ind2;ind3],[ind2;ind3])   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x([ind2;ind3],1),x([ind2;ind3],1),2);          
    %K(ind3,ind3)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind3,1),x(ind3,1),2);              
  elseif i==3
    K = zeros(size(x,1),size(x,1));
    K(ind3,ind3)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind3,1),x(ind3,1),1);      
    %K(ind2,ind2)   = feval(covN1{:}, [location,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),1);              
   elseif i==4
    K = zeros(size(x,1),size(x,1));
    K(ind3,ind3)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind3,1),x(ind3,1),2);  
  elseif i==5
    %K = zeros(size(x,1),size(x,1));
    K = zeros(size(x,1),size(x,1));
    K([ind2;ind3],[ind2;ind3])   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x([ind2;ind3],1),x([ind2;ind3],1),3);
    %K(ind3,ind3)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind3,1),x(ind3,1),3);    
  elseif i==6
    K = zeros(size(x,1),size(x,1));      
    %K(ind1,ind1)   = feval(covN1{:}, [location,steepness1,hyp(4),hyp(5)], x(ind1,1),x(ind1,1),4);
    K([ind2;ind3],[ind2;ind3])   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x([ind2;ind3],1),x([ind2;ind3],1),4);    
    %K(ind3,ind3)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind3,1),x(ind3,1),4);        
  elseif i==7
      %keyboard
      K = zeros(size(x,1),size(x,1));
     K(ind3,ind3)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind3,1),x(ind3,1),3);
  elseif i==8
      K = zeros(size(x,1),size(x,1));
      K(ind3,ind3)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind3,1),x(ind3,1),4);
  elseif i==9
     K = feval('covSEiso',[hyp(9),hyp(10)],x(:,1),x(:,1),1);
  elseif i==10
     K = feval('covSEiso',[hyp(9),hyp(10)],x(:,1),x(:,1),2);
  else
    error('Unknown hyperparameter')
  end
end