function K = covBranchingProcess_4B(hyp, x, z, i)

%A four component branching process (all branch from a latent process)
 
%Note: major scope for speedup. Can predcalculate the distances rather than keep calling.
%Now runs with FITC. 

if nargin<2, K = '26'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists

xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

%Changepoint params.
location1  = hyp(1);  %2
steepness1 = hyp(2);  %2 %Latent 1
location2  = hyp(3);  %2
steepness2 = hyp(4);  %2 %Latent 2
location3  = hyp(5);  %2
steepness3 = hyp(6);  %2 %B1
location4  = hyp(7);  %2
steepness4 = hyp(8);  %2 %B2
location5  = hyp(9);  %2
steepness5 = hyp(10); %2 %B3
location6  = hyp(11); %2
steepness6 = hyp(12); %2 %B4

k1    = {'covMaterniso',3}; %Base kernel
covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}}; %Changepoint kernel

if nargin<4 %Covariances  
    ind1 = find(x(:,2)==1);
    ind2 = find(x(:,2)==2);
    ind3 = find(x(:,2)==3);    
    ind4 = find(x(:,2)==4);        

    if dg    
        K     = feval(k1{:},[hyp(25),hyp(26)],x(:,1)); %Base process
                
       if isempty([ind1;ind2])==0
        K_l1   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)], x([ind1;ind2],1)); %Latent process branch 1/2        
        K([ind1;ind2],[ind1;ind2]) = K([ind1;ind2],[ind1;ind2])+K_l1;
            if isempty(ind1)==0
                K1   = feval(covN1{:},  [location3,steepness3,hyp(17),hyp(18)], x(ind1,1)); %Branch 1
                K(ind1,ind1) = K(ind1,ind1)+K1;
            end        
            if isempty(ind2)==0
                K2   = feval(covN1{:},  [location4,steepness4,hyp(19),hyp(20)], x(ind2,1)); %Branch 2
                K(ind2,ind2) = K(ind2,ind2)+K2;
            end 
       end
        
       
       
        if isempty([ind3;ind4])==0
        K_l2   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)], x([ind3;ind4],1));%Latent process branch 3/4            
        K([ind3;ind4],[ind3;ind4]) = K([ind3;ind4],[ind3;ind4])+K_l2;        
            if isempty(ind3)==0
                K3   = feval(covN1{:},  [location5,steepness5,hyp(21),hyp(22)], x(ind3,1)); %Branch 3
                K(ind3,ind3) = K(ind3,ind3)+K3;
            end        
            if isempty(ind4)==0
                K4   = feval(covN1{:},  [location6,steepness6,hyp(23),hyp(24)], x(ind4,1)); %Branch 4
                K(ind4,ind4) = K(ind4,ind4)+K4;  
            end
        end                                        
        
        K            = diag(K);
    else        
        if  xeqz   
        
         K     = feval(k1{:},[hyp(25),hyp(26)],x(:,1)); %Base process
                
       if isempty([ind1;ind2])==0
        K_l1   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)], x([ind1;ind2],1)); %Latent process branch 1/2        
        K([ind1;ind2],[ind1;ind2]) = K([ind1;ind2],[ind1;ind2])+K_l1;
            if isempty(ind1)==0
                K1   = feval(covN1{:},  [location3,steepness3,hyp(17),hyp(18)], x(ind1,1)); %Branch 1
                K(ind1,ind1) = K(ind1,ind1)+K1;
            end        
            if isempty(ind2)==0
                K2   = feval(covN1{:},  [location4,steepness4,hyp(19),hyp(20)], x(ind2,1)); %Branch 2
                K(ind2,ind2) = K(ind2,ind2)+K2;
            end 
       end
        
        if isempty([ind3;ind4])==0
        K_l2   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)], x([ind3;ind4],1));%Latent process branch 3/4            
        K([ind3;ind4],[ind3;ind4]) = K([ind3;ind4],[ind3;ind4])+K_l2;        
            if isempty(ind3)==0
                K3   = feval(covN1{:},  [location5,steepness5,hyp(21),hyp(22)], x(ind3,1)); %Branch 3
                K(ind3,ind3) = K(ind3,ind3)+K3;
            end        
            if isempty(ind4)==0
                K4   = feval(covN1{:},  [location6,steepness6,hyp(23),hyp(24)], x(ind4,1)); %Branch 4
                K(ind4,ind4) = K(ind4,ind4)+K4;  
            end
        end 
        
        
        
        else %Cross covariances

        ind5 = find(z(:,2)==1);
        ind6 = find(z(:,2)==2);
        ind7 = find(z(:,2)==3);
        ind8 = find(z(:,2)==4);    
         
        K     = feval(k1{:},[hyp(25),hyp(26)],x(:,1),z(:,1)); %Base process
                
       if isempty([ind1;ind2])==0 & isempty([ind5;ind6])==0
        K_l1   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)], x([ind1;ind2],1),z([ind5;ind6],1)); %Latent process branch 1/2        
        K([ind1;ind2],[ind5;ind6]) = K([ind1;ind2],[ind5;ind6])+K_l1;

            if isempty(ind1)==0 & isempty(ind5)==0
                K1   = feval(covN1{:},  [location3,steepness3,hyp(17),hyp(18)], x(ind1,1),z(ind5,1)); %Branch 1
                K(ind1,ind5) = K(ind1,ind5)+K1;
            end        
            if isempty(ind2)==0 & isempty(ind6)==0
                K2   = feval(covN1{:},  [location4,steepness4,hyp(19),hyp(20)], x(ind2,1),z(ind6,1)); %Branch 2
                K(ind2,ind6) = K(ind2,ind6)+K2;
            end 
       end
        
        if isempty([ind3;ind4])==0 & isempty([ind7;ind8])==0
        K_l2   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)], x([ind3;ind4],1),z([ind7;ind8],1));%Latent process branch 3/4            
        K([ind3;ind4],[ind7;ind8]) = K([ind3;ind4],[ind7;ind8])+K_l2;        
            if isempty(ind3)==0 & isempty(ind7)==0
                K3   = feval(covN1{:},  [location5,steepness5,hyp(21),hyp(22)], x(ind3,1),z(ind7,1)); %Branch 3
                K(ind3,ind7) = K(ind3,ind7)+K3;
            end        
            if isempty(ind4)==0 & isempty(ind8)==0
                K4   = feval(covN1{:},  [location6,steepness6,hyp(23),hyp(24)], x(ind4,1),z(ind8,1)); %Branch 4
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

if xeqz==0 & dg==0
    ind5 = find(z(:,2)==1);
    ind6 = find(z(:,2)==2);
    ind7 = find(z(:,2)==3);
    ind8 = find(z(:,2)==4);
end

  if i==1      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        if isempty([ind1;ind2])==0 & isempty([ind5;ind6])==0,K([ind1;ind2],[ind5;ind6])   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)],  x([ind1;ind2],1),z([ind5;ind6],1),1);end
        else
        K = zeros(size(x,1),size(x,1));
        K([ind1;ind2],[ind1;ind2])   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)],  x([ind1;ind2],1),x([ind1;ind2],1),1);      
        if dg,K = diag(K);end
        end

    
  elseif i==2      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind1;ind2])==0 & isempty([ind5;ind6])==0,K([ind1;ind2],[ind5;ind6])   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)],  x([ind1;ind2],1),z([ind5;ind6],1),2);end
        else      
        K = zeros(size(x,1),size(x,1));
        K([ind1;ind2],[ind1;ind2])   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)],  x([ind1;ind2],1),x([ind1;ind2],1),2);      
        if dg,K = diag(K);end
        end
    
 elseif i==3     
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind3;ind4])==0 & isempty([ind7;ind8])==0,K([ind3;ind4],[ind7;ind8])   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)],  x([ind3;ind4],1),z([ind7;ind8],1),1);end
        else
        K = zeros(size(x,1),size(x,1));
        K([ind3;ind4],[ind3;ind4])   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)],  x([ind3;ind4],1),x([ind3;ind4],1),1);      
        if dg,K = diag(K);end
        end
    
  elseif i==4    
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind3;ind4])==0 & isempty([ind7;ind8])==0,K([ind3;ind4],[ind7;ind8])   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)],  x([ind3;ind4],1),z([ind7;ind8],1),2);end
        else
        K = zeros(size(x,1),size(x,1));
        K([ind3;ind4],[ind3;ind4])   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)],  x([ind3;ind4],1),x([ind3;ind4],1),2);              
        if dg,K = diag(K);end
         end        
      
  elseif i==5      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind1])==0 & isempty([ind5])==0,K(ind1,ind5)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind1,1),z(ind5,1),1);end
        else
        K = zeros(size(x,1),size(x,1));
        K(ind1,ind1)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind1,1),x(ind1,1),1);
        if dg,K = diag(K);end
        end
    
  elseif i==6   
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        if isempty([ind1])==0 & isempty([ind5])==0,K(ind1,ind5)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind1,1),z(ind5,1),2);end
        else
        K = zeros(size(x,1),size(x,1));
        K(ind1,ind1)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind1,1),x(ind1,1),2);
        if dg,K = diag(K);end
        end
      
  elseif i==7     
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1)); 
        if isempty([ind2])==0 & isempty([ind6])==0,K(ind2,ind6)   = feval(covN1{:}, [location4,steepness4,hyp(19),hyp(20)], x(ind2,1),z(ind6,1),1);end
        else      
        K = zeros(size(x,1),size(x,1));
        K(ind2,ind2)   = feval(covN1{:}, [location4,steepness4,hyp(19),hyp(20)], x(ind2,1),x(ind2,1),1);
        if dg,K = diag(K);end
        end      

  elseif i==8      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind2])==0 & isempty([ind6])==0,K(ind2,ind6)   = feval(covN1{:}, [location4,steepness4,hyp(19),hyp(20)], x(ind2,1),z(ind6,1),2); end
        else
        K = zeros(size(x,1),size(x,1));
        K(ind2,ind2)   = feval(covN1{:}, [location4,steepness4,hyp(19),hyp(20)], x(ind2,1),x(ind2,1),2);    
        if dg,K = diag(K);end
        end        

elseif i==9    
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind3])==0 & isempty([ind7])==0,K(ind3,ind7)   = feval(covN1{:}, [location5,steepness5,hyp(21),hyp(22)], x(ind3,1),z(ind7,1),1);end
        else
        K = zeros(size(x,1),size(x,1));
        K(ind3,ind3)   = feval(covN1{:}, [location5,steepness5,hyp(21),hyp(22)], x(ind3,1),x(ind3,1),1);
        if dg,K = diag(K);end
        end    

  elseif i==10
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind3])==0 & isempty([ind7])==0,K(ind3,ind7)   = feval(covN1{:}, [location5,steepness5,hyp(21),hyp(22)], x(ind3,1),z(ind7,1),2);end
        else
        K = zeros(size(x,1),size(x,1));
        K(ind3,ind3)   = feval(covN1{:}, [location5,steepness5,hyp(21),hyp(22)], x(ind3,1),x(ind3,1),2);  
        if dg,K = diag(K);end
        end    
    
  elseif i==11  
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind4])==0 & isempty([ind8])==0,K(ind4,ind8)   = feval(covN1{:}, [location6,steepness6,hyp(23),hyp(24)], x(ind4,1),z(ind8,1),1);end
        else
        K = zeros(size(x,1),size(x,1));
        K(ind4,ind4)   = feval(covN1{:}, [location6,steepness6,hyp(23),hyp(24)], x(ind4,1),x(ind4,1),1);
        if dg,K = diag(K);end
        end
      
  elseif i==12
    
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind4])==0 & isempty([ind8])==0,K(ind4,ind8)   = feval(covN1{:}, [location6,steepness6,hyp(23),hyp(24)], x(ind4,1),z(ind8,1),2);end
        else
        K = zeros(size(x,1),size(x,1));
        K(ind4,ind4)   = feval(covN1{:}, [location6,steepness6,hyp(23),hyp(24)], x(ind4,1),x(ind4,1),2);
        if dg,K = diag(K);end
        end      
  
  elseif i==13  
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind1;ind2])==0 & isempty([ind5;ind6])==0,K([ind1;ind2],[ind5;ind6])   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)], x([ind1;ind2],1),z([ind5;ind6],1),3);end
        else
        K = zeros(size(x,1),size(x,1));
        K([ind1;ind2],[ind1;ind2])   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)], x([ind1;ind2],1),x([ind1;ind2],1),3);      
        if dg,K = diag(K);end
        end      
  
  elseif i==14   
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind1;ind2])==0 & isempty([ind5;ind6])==0,K([ind1;ind2],[ind5;ind6])   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)], x([ind1;ind2],1),z([ind5;ind6],1),4);end 
        else
        K = zeros(size(x,1),size(x,1));
        K([ind1;ind2],[ind1;ind2])   = feval(covN1{:}, [location1,steepness1,hyp(13),hyp(14)], x([ind1;ind2],1),x([ind1;ind2],1),4);             
        if dg,K = diag(K);end
        end
      
  elseif i==15    
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind3;ind4])==0 & isempty([ind7;ind8])==0,K([ind3;ind4],[ind7;ind8])   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)], x([ind3;ind4],1),z([ind7;ind8],1),3);end
        else
        K = zeros(size(x,1),size(x,1));
        K([ind3;ind4],[ind3;ind4])   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)], x([ind3;ind4],1),x([ind3;ind4],1),3);                
        if dg,K = diag(K);end
        end      

  elseif i==16      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind3;ind4])==0 & isempty([ind7;ind8])==0,K([ind3;ind4],[ind7;ind8])   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)], x([ind3;ind4],1),z([ind7;ind8],1),4);end
        else      
        K = zeros(size(x,1),size(x,1));
        K([ind3;ind4],[ind3;ind4])   = feval(covN1{:}, [location2,steepness2,hyp(15),hyp(16)], x([ind3;ind4],1),x([ind3;ind4],1),4);        
        if dg,K = diag(K);end
        end
      
  elseif i==17      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind1])==0 & isempty([ind5])==0,K(ind1,ind5)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind1,1),z(ind5,1),3);end
        else
        K = zeros(size(x,1),size(x,1));
        K(ind1,ind1)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind1,1),x(ind1,1),3);
        if dg,K = diag(K);end
        end

  elseif i==18      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind1])==0 & isempty([ind5])==0,K(ind1,ind5)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind1,1),z(ind5,1),4);end
        else
        K = zeros(size(x,1),size(x,1));      
        K(ind1,ind1)   = feval(covN1{:}, [location3,steepness3,hyp(17),hyp(18)], x(ind1,1),x(ind1,1),4);
        if dg,K = diag(K);end
        end

  elseif i==19      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind2])==0 & isempty([ind6])==0,K(ind2,ind6)   = feval(covN1{:}, [location4,steepness4,hyp(19),hyp(20)], x(ind2,1),z(ind6,1),3);end
        else
        K = zeros(size(x,1),size(x,1));
        K(ind2,ind2)   = feval(covN1{:}, [location4,steepness4,hyp(19),hyp(20)], x(ind2,1),x(ind2,1),3);
        if dg,K = diag(K);end
        end

  elseif i==20      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind2])==0 & isempty([ind6])==0,K(ind2,ind6)   = feval(covN1{:}, [location4,steepness4,hyp(19),hyp(20)], x(ind2,1),z(ind6,1),4);end
        else
        K = zeros(size(x,1),size(x,1));      
        K(ind2,ind2)   = feval(covN1{:}, [location4,steepness4,hyp(19),hyp(20)], x(ind2,1),x(ind2,1),4);
        if dg,K = diag(K);end
        end

  elseif i==21
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind3])==0 & isempty([ind7])==0,K(ind3,ind7)   = feval(covN1{:}, [location5,steepness5,hyp(21),hyp(22)], x(ind3,1),z(ind7,1),3);end
        else
        K = zeros(size(x,1),size(x,1));
        K(ind3,ind3)   = feval(covN1{:}, [location5,steepness5,hyp(21),hyp(22)], x(ind3,1),x(ind3,1),3);
        if dg,K = diag(K);end
        end

  elseif i==22    
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind3])==0 & isempty([ind7])==0,K(ind3,ind7)   = feval(covN1{:}, [location5,steepness5,hyp(21),hyp(22)], x(ind3,1),z(ind7,1),4);end
        else
        K = zeros(size(x,1),size(x,1));      
        K(ind3,ind3)   = feval(covN1{:}, [location5,steepness5,hyp(21),hyp(22)], x(ind3,1),x(ind3,1),4);
        if dg,K = diag(K);end
        end

  elseif i==23
    
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind4])==0 & isempty([ind8])==0,K(ind4,ind8)   = feval(covN1{:}, [location6,steepness6,hyp(23),hyp(24)], x(ind4,1),z(ind8,1),3);end
        else
        K = zeros(size(x,1),size(x,1));      
        K(ind4,ind4)   = feval(covN1{:}, [location6,steepness6,hyp(23),hyp(24)], x(ind4,1),x(ind4,1),3);      
        if dg,K = diag(K);end
        end

  elseif i==24      
        if xeqz==0 & dg==0
        K = zeros(size(x,1),size(z,1));
        if isempty([ind4])==0 & isempty([ind8])==0,K(ind4,ind8)   = feval(covN1{:}, [location6,steepness6,hyp(23),hyp(24)], x(ind4,1),z(ind8,1),4);end
        else
        K = zeros(size(x,1),size(x,1));      
        K(ind4,ind4)   = feval(covN1{:}, [location6,steepness6,hyp(23),hyp(24)], x(ind4,1),x(ind4,1),4);          
        if dg,K = diag(K);end
        end

  elseif i==25     
        if xeqz==0 & dg==0
         K = feval(k1{:},[hyp(25),hyp(26)],x(:,1),z(:,1),1);
        else
        K = feval(k1{:},[hyp(25),hyp(26)],x(:,1),x(:,1),1); %Check
        if dg,K = diag(K);end
        end

  elseif i==26
        if xeqz==0 & dg==0
        K = feval(k1{:},[hyp(25),hyp(26)],x(:,1),z(:,1),2);
        else
        K = feval(k1{:},[hyp(25),hyp(26)],x(:,1),x(:,1),2); %Check
        if dg,K = diag(K);end
        end
  else
    error('Unknown hyperparameter')
  end
end