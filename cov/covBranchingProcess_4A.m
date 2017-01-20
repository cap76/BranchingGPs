function K = covBranchingProcess_4A(hyp, x, z, i)

%A four component branching process (all branch from a latent process). Does not 
%yet work with FITC.

if nargin<2, K = '20'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists

xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

%Changepoint params.
location1  = hyp(1);
steepness1 = hyp(2);
location2  = hyp(3);
steepness2 = hyp(4);
location3  = hyp(5);
steepness3 = hyp(6);

k1    = {'covMaterniso',3}; %Base kernel
covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}}; %Changepoint kernel


if nargin<4 %Covariances  
    ind1 = find(x(:,2)==1);
    ind2 = find(x(:,2)==2);
    ind3 = find(x(:,2)==3);    
    ind4 = find(x(:,2)==4);        

    if dg    
        K     = feval(k1{:},[hyp(19),hyp(20)],x(:,1)); %Base process
                
       if isempty([ind1;ind2])==0
        K_l1   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x([ind1;ind2],1)); %Latent process branch 1/2        
        K([ind1;ind2],[ind1;ind2]) = K([ind1;ind2],[ind1;ind2])+K_l1;
        
            if isempty(ind1)==0
                K1   = feval(covN1{:},  [location2,steepness2,hyp(11),hyp(12)], x(ind1,1)); %Branch 1
                K(ind1,ind1) = K(ind1,ind1)+K1;
            end        
            
            if isempty(ind2)==0
                K2   = feval(covN1{:},  [location2,steepness2,hyp(13),hyp(14)], x(ind2,1)); %Branch 2
                K(ind2,ind2) = K(ind2,ind2)+K2;
            end 
       end
        
        if isempty([ind3;ind4])==0
        K_l2   = feval(covN1{:}, [location1,steepness1,hyp(9),hyp(10)], x([ind3;ind4],1));%Latent process branch 3/4            
        K([ind3;ind4],[ind3;ind4]) = K([ind3;ind4],[ind3;ind4])+K_l2;        
            
            if isempty(ind3)==0
                K3   = feval(covN1{:},  [location3,steepness3,hyp(15),hyp(16)], x(ind3,1)); %Branch 3
                K(ind3,ind3) = K(ind3,ind3)+K3;
            end        
            if isempty(ind4)==0
                K4   = feval(covN1{:},  [location3,steepness3,hyp(17),hyp(18)], x(ind4,1)); %Branch 4
                K(ind4,ind4) = K(ind4,ind4)+K4;  
            end
        end                                        
        
        K            = diag(K);
    else        
        if  xeqz   
        
         K     = feval(k1{:},[hyp(19),hyp(20)],x(:,1)); %Base process
                
       if isempty([ind1;ind2])==0
        K_l1   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x([ind1;ind2],1)); %Latent process branch 1/2        
        K([ind1;ind2],[ind1;ind2]) = K([ind1;ind2],[ind1;ind2])+K_l1;
            if isempty(ind1)==0
                K1   = feval(covN1{:},  [location2,steepness2,hyp(11),hyp(12)], x(ind1,1)); %Branch 1
                K(ind1,ind1) = K(ind1,ind1)+K1;
            end        
            if isempty(ind2)==0
                K2   = feval(covN1{:},  [location2,steepness2,hyp(13),hyp(14)], x(ind2,1)); %Branch 2
                K(ind2,ind2) = K(ind2,ind2)+K2;
            end 
       end
        
        if isempty([ind3;ind4])==0
        K_l2   = feval(covN1{:}, [location1,steepness1,hyp(9),hyp(10)], x([ind3;ind4],1));%Latent process branch 3/4            
        K([ind3;ind4],[ind3;ind4]) = K([ind3;ind4],[ind3;ind4])+K_l2;        
            if isempty(ind3)==0
                K3   = feval(covN1{:},  [location3,steepness3,hyp(15),hyp(16)], x(ind3,1)); %Branch 3
                K(ind3,ind3) = K(ind3,ind3)+K3;
            end        
            if isempty(ind4)==0
                K4   = feval(covN1{:},  [location3,steepness3,hyp(17),hyp(18)], x(ind4,1)); %Branch 4
                K(ind4,ind4) = K(ind4,ind4)+K4;  
            end
        end 
        
        else %Cross covariances
        ind5 = find(z(:,2)==1);
        ind6 = find(z(:,2)==2);
        ind7 = find(z(:,2)==3);
        ind8 = find(z(:,2)==4);                
        
        K     = feval(k1{:},[hyp(19),hyp(20)],x(:,1),z(:,1)); %Base process
                
       if isempty([ind1;ind2])==0 & isempty([ind5;ind6])==0
        K_l1   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x([ind1;ind2],1),z([ind5;ind6],1)); %Latent process branch 1/2        
        K([ind1;ind2],[ind5;ind6]) = K([ind1;ind2],[ind5;ind6])+K_l1;

            if isempty(ind1)==0 & isempty(ind5)==0
                K1   = feval(covN1{:},  [location2,steepness2,hyp(11),hyp(12)], x(ind1,1),z(ind5,1)); %Branch 1
                K(ind1,ind5) = K(ind1,ind5)+K1;
            end        
            if isempty(ind2)==0 & isempty(ind6)==0
                K2   = feval(covN1{:},  [location2,steepness2,hyp(13),hyp(14)], x(ind2,1),z(ind6,1)); %Branch 2
                K(ind2,ind6) = K(ind2,ind6)+K2;
            end 
       end
        
        if isempty([ind3;ind4])==0 & isempty([ind7;ind8])==0
        K_l2   = feval(covN1{:}, [location1,steepness1,hyp(9),hyp(10)], x([ind3;ind4],1),z([ind7;ind8],1));%Latent process branch 3/4            
        K([ind3;ind4],[ind7;ind8]) = K([ind3;ind4],[ind7;ind8])+K_l2;        
            if isempty(ind3)==0 & isempty(ind7)==0
                K3   = feval(covN1{:},  [location3,steepness3,hyp(15),hyp(16)], x(ind3,1),z(ind7,1)); %Branch 3
                K(ind3,ind7) = K(ind3,ind7)+K3;
            end        
            if isempty(ind4)==0 & isempty(ind8)==0
                K4   = feval(covN1{:},  [location3,steepness3,hyp(17),hyp(18)], x(ind4,1),z(ind8,1)); %Branch 4
                K(ind4,ind8) = K(ind4,ind8)+K4;  
            end
        end                        
              
        end
end
    
else %Get derivatives (turn this into a loop for arbitrary branching)?
ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);
ind3 = find(x(:,2)==3);
ind4 = find(x(:,2)==4);

  if i==1
    K = zeros(size(x,1),size(x,1));
    K([ind1;ind2],[ind1;ind2])   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)],  x([ind1;ind2],1),x([ind1;ind2],1),1);      
    K([ind3;ind4],[ind3;ind4])   = feval(covN1{:}, [location1,steepness1,hyp(9),hyp(10)], x([ind3;ind4],1),x([ind3;ind4],1),1);      
  elseif i==2
    K = zeros(size(x,1),size(x,1));
    K([ind1;ind2],[ind1;ind2])   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)],  x([ind1;ind2],1),x([ind1;ind2],1),2);      
    K([ind3;ind4],[ind3;ind4])   = feval(covN1{:}, [location1,steepness1,hyp(9),hyp(10)], x([ind3;ind4],1),x([ind3;ind4],1),2);     
  elseif i==3      
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location2,steepness2,hyp(11),hyp(12)], x(ind1,1),x(ind1,1),1);
    K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(13),hyp(14)], x(ind2,1),x(ind2,1),1);    
   elseif i==4
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location2,steepness2,hyp(11),hyp(12)], x(ind1,1),x(ind1,1),2);
    K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(13),hyp(14)], x(ind2,1),x(ind2,1),2);     
  elseif i==5
    K = zeros(size(x,1),size(x,1));
    K(ind3,ind3)   = feval(covN1{:}, [location3,steepness3,hyp(15),hyp(16)], x(ind3,1),x(ind3,1),1);
    K(ind4,ind4)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind4,1),x(ind4,1),1);    
   elseif i==6
    K = zeros(size(x,1),size(x,1));
    K(ind3,ind3)   = feval(covN1{:}, [location3,steepness3,hyp(15),hyp(16)], x(ind3,1),x(ind3,1),2);
    K(ind4,ind4)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind4,1),x(ind4,1),2);         
  elseif i==7
    K = zeros(size(x,1),size(x,1));
    K([ind1;ind2],[ind1;ind2])   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x([ind1;ind2],1),x([ind1;ind2],1),3);      
   elseif i==8
    K = zeros(size(x,1),size(x,1));
    K([ind1;ind2],[ind1;ind2])   = feval(covN1{:}, [location1,steepness1,hyp(7),hyp(8)], x([ind1;ind2],1),x([ind1;ind2],1),4);             
   elseif i==9 
    K = zeros(size(x,1),size(x,1));
    K([ind3;ind4],[ind3;ind4])   = feval(covN1{:}, [location1,steepness1,hyp(9),hyp(10)], x([ind3;ind4],1),x([ind3;ind4],1),3);                
   elseif i==10
    K = zeros(size(x,1),size(x,1));
    K([ind3;ind4],[ind3;ind4])   = feval(covN1{:}, [location1,steepness1,hyp(9),hyp(10)], x([ind3;ind4],1),x([ind3;ind4],1),4);        
  elseif i==11
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location2,steepness2,hyp(11),hyp(12)], x(ind1,1),x(ind1,1),3);
  elseif i==12
    K = zeros(size(x,1),size(x,1));      
    K(ind1,ind1)   = feval(covN1{:}, [location2,steepness2,hyp(11),hyp(12)], x(ind1,1),x(ind1,1),4);
  elseif i==13
    K = zeros(size(x,1),size(x,1));
  K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(13),hyp(14)], x(ind2,1),x(ind2,1),3);
  elseif i==14
    K = zeros(size(x,1),size(x,1));      
  K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(13),hyp(14)], x(ind2,1),x(ind2,1),4);
  elseif i==15
    K = zeros(size(x,1),size(x,1));
    K(ind3,ind3)   = feval(covN1{:}, [location3,steepness3,hyp(15),hyp(16)], x(ind3,1),x(ind3,1),3);
  elseif i==16
    K = zeros(size(x,1),size(x,1));      
    K(ind3,ind3)   = feval(covN1{:}, [location3,steepness3,hyp(15),hyp(16)], x(ind3,1),x(ind3,1),4);
 elseif i==17
    K = zeros(size(x,1),size(x,1));      
    K(ind4,ind4)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind4,1),x(ind4,1),3);      
 elseif i==18
    K = zeros(size(x,1),size(x,1));      
    K(ind4,ind4)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind4,1),x(ind4,1),4);          
  elseif i==19
     K = feval(k1{:},[hyp(19),hyp(20)],x(:,1),x(:,1),1); %Check
  elseif i==20
     K = feval(k1{:},[hyp(19),hyp(20)],x(:,1),x(:,1),2); %Check
  else
    error('Unknown hyperparameter')
  end
end