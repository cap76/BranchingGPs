function K = covBranchingProcess_5Rec(hyp, x, z, i)

%A simple branching process, in which a process Y(t) branches from the underlying process
%X(t) at some time. Prior to this P(X(t))=P(Y(t)).

%Note: should now fully function within gpml including FITC. 

if nargin<2, K = '16'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

location1   = hyp(1); %Branch point1
steepness1 = hyp(2); %Smoothness in branching
location2   = hyp(3); %Rec point1
steepness2 = hyp(4); %Smoothness in recomb.
location3   = hyp(5); %Branch point2
steepness3 = hyp(6); %Smoothness in branching2
location4   = hyp(7); %Branch point
steepness4 = hyp(8); %Smoothness in branching

k1={'covMaterniso',3}; %Base kernel


covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}}; %Branch process
covN2 = {@covChangePointMultiD, {1, {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}}, @covZero}}; %Rec process


if nargin<4                                                        % covariances  
    ind1 = find(x(:,2)==1); %Branc rec
    ind2 = find(x(:,2)==2); %2nd branch index
    ind3 = find(x(:,2)==3); %2nd branch index    
    
     if dg    
         
        K    = feval(k1{:},[hyp(15),hyp(16)],x(:,1)); %Base branch
        
        if isempty(ind1)==0 
        K1   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1));
        K(ind1,ind1) = K(ind1,ind1) + K1;        
        end
        
        if isempty(ind2)==0 
        K2   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1));
        K(ind2,ind2) = K(ind2,ind2) + K2;        
        end        
        
        if isempty(ind3)==0 
        K3   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1));
        K(ind3,ind3) = K(ind3,ind3) + K3;        
        end
        
        K    = diag(K);        %Faster way of computing diag. 

     else        
        if  xeqz   %Symmetric
        K    = feval(k1{:},[hyp(15),hyp(16)],x(:,1)); %Faster way of computing this if on a grid.
            
        if isempty(ind1)==0 
        K1   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1));
        K(ind1,ind1) = K(ind1,ind1) + K1;        
        end
        
        if isempty(ind2)==0 
        K2   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1));
        K(ind2,ind2) = K(ind2,ind2) + K2;        
        end        
        
        if isempty(ind3)==0 
        K3   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1));
        K(ind3,ind3) = K(ind3,ind3) + K3;        
        end
            
        else %Cross covariances
            ind4 = find(z(:,2)==1);
            ind5 = find(z(:,2)==2);        
            ind6 = find(z(:,2)==3);        

            K    = feval(k1{:},[hyp(15),hyp(16)],x(:,1),z(:,1));
            
            if isempty(ind1)==0 & isempty(ind4)==0
            K1   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1), z(ind4,1));
            K(ind1,ind4) = K(ind1,ind4) + K1;        
            end
        
            if isempty(ind2)==0 & isempty(ind5)==0
            K2   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1), z(ind5,1));
            K(ind2,ind5) = K(ind2,ind5) + K2;        
            end        
        
            if isempty(ind3)==0 & isempty(ind6)==0
            K3   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1), z(ind6,1));
            K(ind3,ind6) = K(ind3,ind6) + K3;        
            end           
            
            
        end    
    end
    
else
    
ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);
ind3 = find(x(:,2)==3);
if xeqz==0 & dg==0
    ind4 = find(z(:,2)==1);
    ind5 = find(z(:,2)==2);
    ind6 = find(z(:,2)==3);
end

  % derivatives       
 if i==1
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind1)==0 & isempty(ind4)==0,K(ind1,ind4)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1),z(ind4,1),1);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind1)==0,K(ind1,ind1)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1),x(ind1,1),1);end      
        if dg,K = diag(K);end
        end                
        
 elseif i==2
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind1)==0 & isempty(ind4)==0,K(ind1,ind4)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1),z(ind4,1),2);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind1)==0,K(ind1,ind1)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1),x(ind1,1),2);end      
        if dg,K = diag(K);end
        end         
        
 elseif i==3
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        if isempty(ind2)==0 & isempty(ind5)==0,K(ind2,ind5)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind2,1),z(ind5,1),3);end 
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind2)==0,K(ind2,ind2)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind2,1),x(ind2,1),3);end          
        if dg,K = diag(K);end
        end
        
 elseif i==4
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind1)==0 & isempty(ind4)==0,K(ind1,ind4)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1),z(ind4,1),4);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind1)==0,K(ind1,ind1)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1),x(ind1,1),4);end      
        if dg,K = diag(K);end
        end    
        
 elseif i==5
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind2)==0 & isempty(ind5)==0,K(ind2,ind5)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1),z(ind5,1),1);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1),x(ind2,1),1);end      
        if dg,K = diag(K);end
        end            
         
 elseif i==6
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind2)==0 & isempty(ind5)==0,K(ind2,ind5)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1),z(ind5,1),2);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1),x(ind2,1),2);end      
        if dg,K = diag(K);end
        end            
                        
 elseif i==7
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind3)==0 & isempty(ind6)==0,K(ind3,ind6)   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1),z(ind6,1),1);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind3)==0,K(ind3,ind3)   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1),x(ind3,1),1);end      
        if dg,K = diag(K);end
        end            
 
        
 elseif i==8
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind3)==0 & isempty(ind6)==0,K(ind3,ind6)   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1),z(ind6,1),2);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind3)==0,K(ind3,ind3)   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1),x(ind3,1),2);end      
        if dg,K = diag(K);end
        end            
                                
 elseif i==9
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind1)==0 & isempty(ind4)==0,K(ind1,ind4)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1),z(ind4,1),5);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind1)==0,K(ind1,ind1)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1),x(ind1,1),5);end      
        if dg,K = diag(K);end
        end       
        
 elseif i==10
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind1)==0 & isempty(ind4)==0,K(ind1,ind4)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1),z(ind4,1),6);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind1)==0,K(ind1,ind1)   = feval(covN2{:}, [location1,steepness1,location2,steepness2,hyp(9),hyp(10)], x(ind1,1),x(ind1,1),6);end      
        if dg,K = diag(K);end
        end           
                
 elseif i==11
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind2)==0 & isempty(ind5)==0,K(ind2,ind5)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1),z(ind5,1),3);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1),x(ind2,1),3);end      
        if dg,K = diag(K);end
        end                     

elseif i==12
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind2)==0 & isempty(ind5)==0,K(ind2,ind5)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1),z(ind5,1),4);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind2)==0,K(ind2,ind2)   = feval(covN1{:}, [location3,steepness3,hyp(11),hyp(12)], x(ind2,1),x(ind2,1),4);end      
        if dg,K = diag(K);end
        end           
        
 elseif i==13
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind3)==0 & isempty(ind6)==0,K(ind3,ind6)   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1),z(ind6,1),3);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind3)==0,K(ind3,ind3)   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1),x(ind3,1),3);end      
        if dg,K = diag(K);end
        end            

        
 elseif i==14
      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        %If indexes exist
        if isempty(ind3)==0 & isempty(ind6)==0,K(ind3,ind6)   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1),z(ind6,1),4);end
        else
        K = zeros(size(x,1),size(x,1));
        if isempty(ind3)==0,K(ind3,ind3)   = feval(covN1{:}, [location4,steepness4,hyp(13),hyp(14)], x(ind3,1),x(ind3,1),4);end      
        if dg,K = diag(K);end
        end            
         
        
    elseif i==15
        
        if xeqz==0 & dg==0
        K = feval(k1{:},[hyp(15),hyp(16)],x(:,1),z(:,1),1);
        else
        K = feval(k1{:},[hyp(15),hyp(16)],x(:,1),x(:,1),1);
        if dg,K = diag(K);end
        end
        
     elseif i==16
        if xeqz==0 & dg==0
        K = feval(k1{:},[hyp(15),hyp(16)],x(:,1),z(:,1),2);
        else
        K = feval(k1{:},[hyp(15),hyp(16)],x(:,1),x(:,1),2);
        if dg,K = diag(K);end
        end
        
 else      
    error('Unknown hyperparameter')
  end
end