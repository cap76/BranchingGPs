function K = covBranchingProcess_2A(hyp, x, z, i)

%A simple branching process, in which a process Y(t) branches from the underlying process
%X(t) at some time. Prior to this P(X(t))=P(Y(t)).

if nargin<2, K = '6'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

location   = hyp(1); %Branch point
steepness1 = hyp(2); %Smoothness in branching
covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
  
if nargin<4                                                        % covariances  
    ind1 = find(x(:,2)==1); %Main branch index
    ind2 = find(x(:,2)==2); %2nd branch index
    
     if dg    
        K    = covSEiso([hyp(5),hyp(6)],x(:,1)); %Base branch
        K2   = feval(covN1{:}, [location,steepness1,hyp(3),hyp(4)], x(ind2,1));
        K(ind2,ind2) = K(ind2,ind2) + K2;        
        K    = diag(K);
     else        
    if  xeqz   
    K    = covSEiso([hyp(5),hyp(6)],x(:,1));
    K2   = feval(covN1{:}, [location,steepness1,hyp(3),hyp(4)], x(ind2,1));
    K(ind2,ind2) = K(ind2,ind2) + K2;
    else %Cross covariances
        ind3 = find(z(:,2)==1);
        ind4 = find(z(:,2)==2);        
        K    = covSEiso([hyp(5),hyp(6)],x(:,1),z(:,1));
        K2   = feval(covN1{:}, [location,steepness1,hyp(3),hyp(4)], x(ind2,1),z(ind4,1));
        if isempty(ind4)==0
        K(ind2,ind4) = K(ind2,ind4)+K2;
        end
    end    
end
    
else
    
ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);
  % derivatives
  if i==1
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),1);      
  elseif i==2
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),2);          
  elseif i==3
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),3);
  elseif i==4
    K = zeros(size(x,1),size(x,1));      
    K(ind2,ind2)   = feval(covN1{:}, [location,steepness1,hyp(3),hyp(4)], x(ind2,1),x(ind2,1),4);    
  elseif i==5
     K = feval('covSEiso',[hyp(5),hyp(6)],x(:,1),x(:,1),1);
  elseif i==6
     K = feval('covSEiso',[hyp(5),hyp(6)],x(:,1),x(:,1),2);
  else
    error('Unknown hyperparameter')
  end
end