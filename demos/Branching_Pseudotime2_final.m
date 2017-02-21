function Output = Branching_Pseudotime2_final(seed,missingT);

%Add path to various toolboxes (gpml, netlab and changepoint kernels)
addpath(genpath('../'))
%addpath(genpath('../Pseudotime_Algorithms/Functions/Others'))
%addpath(genpath('../Pseudotime_Algorithms/Functions/gpml-matlab-v3.6-2015-07-07/'))
%addpath(genpath('../Code/gpss-research/source/gpml/cov/'))
%addpath(genpath('../Code/gpss-research/source/matlab/custom-cov/'))
%addpath(genpath('../Code/deepGP/netlab3_3/'))

%Example covariance functions covBranching_v5a, covBranching_v5b,
%covBranching_v4b. 2 branch and 3 branch processes.

try
    sdadasd load(['/Users/christopher_penfold/Desktop/Branching_GPs/PseudotimeCluster/New_new_new/Marker_Pseudotime_' num2str(seed) '_' num2str(missingT) '_2E.mat'])
    %load(['Marker_Pseudotime_' num2str(seed) '_' num2str(missingT) '_2E.mat'])
    rng(Output.s)
    if isempty(Output.stepno)
        Output.stepno = size(Output.L,1);
    end
    StartNo = Output.stepno+1;
    NoMCMC  = 20000;
    Data    = Output.Data;
    L       = Output.L;
    Prt     = Output.Pr;
    Param   = Output.Param;
    Xstar   = Output.Xstar;
    count   = Data.update.count;
    uc1 = Data.update.uc1;
    uc2 = Data.update.uc2;
    uc3 = Data.update.uc3;
    uc4 = Data.update.uc4;
    uc5 = Data.update.uc5;
    uc6 = Data.update.uc6;
    uc7 = Data.update.uc7;
    ub1 = Data.update.ub1;
    ub4 = Data.update.ub4;    
    t1  = Data.orig.t';
    options = Data.orig.options;
catch

%Reset random number generator
rng(seed)
NoMCMC = 20000; %No. steps in MCMC

%Grab some data ...
[D1,t1,Type,Sex,uID,FixLabel,FixTime,Assign,nocells] = initialisedata(missingT);
%Initialise the prior distributions for all the hyperaprams
[options, prior, altprior] = initialiseprior(D1); %Initialise prior distributions
%Initialise the modes (assign unlabelled data on the fly)


[Data,X,Y,Xstar] = initialiseBranch(D1.data(4:end,:),t1,Type,Sex,uID,FixLabel,FixTime,Assign,prior);

%Tune switch points by fix other hyperparameters from here on out.
%Data.orig.options = options;
Data.update.im{3} = altprior;

Data.update.AR1 = [];
Data.update.AR2 = [];

%Initilise some storage
L   = zeros(NoMCMC,size(Data.update.y,2)); %Likelihoods
Prt = zeros(NoMCMC,size(Data.update.y,2)); %Prior over times

%Get initial likelihoods
for arc = 1:size(Data.update.y,2)              %Inference method
    ind = find((Data.update.y(:,arc))~=0);    %In each GP model ignore missing data ~=0 (change to pi when want to include everything).
    par = {'meanConst','covBranchingProcess_2E','likGauss',Data.update.x(ind,:),Data.update.y(ind,arc)};
    [L(1,arc) dL] = gp(Data.update.hyp{arc}, Data.update.im, par{:});
end
%keyboard
%Log prior probability over timepoints
for arc= 1:length(Data.update.origt)    
    if  Data.update.origt(arc)==min(Data.update.origt) %Uniform prior for ESC/anything with no capture time
        Prt(1,arc) = feval(Data.update.p1{:},Data.update.x(arc,1));
    else                                             %Everything else is Gaussian centred on capture time
        pri = {@priorGauss,Data.update.origt(arc),.01};
        Prt(1,arc) = feval(pri{:},Data.update.x(arc,1));
    end    
end

%Param.Store = cell(NoMCMC,1);
Param.Store{1} = Data; %Store everything at each step. Can now vectorise now we've sorted out a bunch of covariance functions.
options(15)    = 200;  %Length of HMC for each hyperparameter updates
count          = 1;    %Acceptence count (for tuning some hyperparams)
uc1 = 0; uc2 = 0; uc3 = 0;  uc4 = 0; uc5 = 0; uc6 = 0; uc7 = 0; ub1 = 0; ub4 = 0; %Other counts 

Data.orig.options = options;
Data.update.uc1 = uc1;
Data.update.uc2 = uc2;
Data.update.uc3 = uc3;
Data.update.uc4 = uc4;
Data.update.uc5 = uc5;
Data.update.uc6 = uc6;
Data.update.uc7 = uc7;
Data.update.ub1 = ub1;
Data.update.ub4 = ub4;

Data.update.counts = [uc1,uc2,uc3,uc4,uc5,uc6,uc7];
StartNo = 2;

end


for i = StartNo:NoMCMC   
    
 a = rand(1,1);
    
    if a<0.5    
    %Gaussian perturbation of the pseudotime time for a subset of cells    
        %if rand(1,1)<0.5 %Update many genes by small amounts
            [Data,L(i,:),Prt(i,:),u1] = updateTimeGauss(Data,L(i-1,:),Prt(i-1,:));    
            uc1 = uc1 + abs(u1); %Count accepted
            ub1 = ub1+1;         %Count tries
        %else
        %    %Global cell updade. A few cells bigger jumps.
        %    [Data,L(i,:),Prt(i,:),u2] = updateGlobalTimeGauss(Data,L(i-1,:),Prt(i-1,:));
        %    uc2 = uc2+abs(u2);
        %end    
    elseif a>=0.5 & a<0.75
        %if rand(1,1)<0.5
            %Swap the pseudotime of adjacent cells (in a given branch)
            [Data,L(i,:),Prt(i,:),u3] = updateTimeSwap(Data,L(i-1,:),Prt(i-1,:));
            uc3 = uc3 + abs(u3);       
        %else
        %    %Globally swap cells (in a branch). Maybe between branches?
        %    [Data,L(i,:),Prt(i,:),u4] = updateGlobalTimeSwap(Data,L(i-1,:),Prt(i-1,:));    
        %    uc4 = uc4+abs(u4);
        %    ub4 = ub4+1; 
        %end
    else%if a>=0.5 & a<0.75        
        Prt(i,:) = Prt(i-1,:);
        %if rand(1,1)<0.5
            %Update the branch assignment of a single cell
            [Data,L(i,:),u5] = updateAssignment(Data,L(i-1,:));
            uc5 = uc5 + abs(u5);    
        %else
        %    %Try flipping a random set of branch assignments
        %    [Data,L(i,:),u6] = updateGlobalAssignment(Data,L(i-1,:));
        %    uc6 = uc6 + abs(u6);                           
        %end        
    %else        
        %if rand(1,1)<0.5
     %   Prt(i,:) = Prt(i-1,:);
     %   [Data,L(i,:),u7] = updateOutliers(Data,L(i-1,:),Prt(i,:));        
        %else
        %Marginally assign bunch.
        %[Data,L(i,:),u8] = marginalupdateOutliers(Data,L(i-1,:));        
        %uc6 = uc* + abs(u8);
        %ub8 = ub8+1; 
        %end
      %  
    end    
        
    %Tune some parameters every 100 steps.
    if double(int64(i/100))==(i/100)
        AR = uc1./(ub1); %Acceptence rate of the last 100 samples
        
        if abs(AR-0.25)<=0.1 %Keep the same
            
        elseif abs(AR-0.25)>0.1 & AR>0.25   %Accepting too many (increase number to swap/variance of step)
            Data.update.TR = Data.update.TR*1.1;
            Data.update.p  = Data.update.p*1.1;                       
            Data.update.TR = min(Data.update.TR,1);
            Data.update.p  = min(Data.update.p,1);
        else             %Accepting too few (decrease number and variance of step)
            Data.update.TR = Data.update.TR*0.9;
            Data.update.p  = Data.update.p*0.9;     
            Data.update.TR = max(Data.update.TR,0.001);
            Data.update.p  = max(Data.update.p,0.001);
        end           
         
        Data.update.AR1 = [Data.update.AR1;AR];
        
         AR = uc4./(ub4); %Acceptence rate of the last 100 samples
        
        if abs(AR-0.25)<=0.1 %Keep the same
            
        elseif abs(AR-0.25)>0.1 & AR>0.25   %Accepting too many (increase number to swap/variance of step)
            Data.update.pgs = Data.update.pgs*1.1;
            Data.update.pgs  = min(Data.update.pgs,1);
        else             %Accepting too few (decrease number and variance of step)
            Data.update.pgs  = Data.update.pgs*0.9;                 
            Data.update.pgs  = max(Data.update.pgs,0.001);            
        end           
        
        Data.update.counts = Data.update.counts+[uc1,uc2,uc3,uc4,uc5,uc6,uc7];        
        Data.update.AR2 = [Data.update.AR2;AR];
        uc1 = 0; uc2 = 0; uc3 = 0;  uc4 = 0; uc5 = 0; uc6 = 0; uc7 = 0; ub1 = 0; ub4 = 0; %Other counts         
    end           
    
     %if i==Data.update.burnin
     %Oooh ... after burn-in lets update all labels! Reset the "fixed labels"
     %   FixLabel = zeros(1,length(t1));
     %   Data.update.origFixLabel = FixLabel;        
     %   Data.update.FixLabel = FixLabel;
     %   
     %   %Wait ... don't update the bulk data ...
     %   %Data.update.FixLabel(1,find(Data.update.uID)>nocells) = ones;               
     %end
    
     %Any other types of update? Perhaps think about more global
     %(heuristic updates with a MH type acceptentence) i.e., branch
     %switch. Complete switching of branch labels?
       
     %Every 100 steps do HMC on hyperparams. Maybe break this down ... the branch times are what we want to update really. 
         
         
         %Next test the optimisation of branch times/rates only. All other
         %hyperparms fixed.
         if double(int64(i/100))==(i/100) & i>Data.update.burnin
             
         count = count+1;
         for arc = 1:size(Data.update.y,2)
         ind = find((Data.update.y(:,arc)~=0));              
         par = {'meanConst','covBranchingProcess_2E','likGauss',Data.update.x(ind,:),Data.update.y(ind,arc),Data.update.hyp{arc}};
         [samples, energies] = hmc2('gpfunc', unwrap(Data.update.hyp{arc}), options, 'gpgrad', Data.update.im, par{:});     
         hypupdate  = rewrap(Data.update.hyp{arc},samples(end,:));
         Data.update.hyp{arc} = rewrap(Data.update.hyp{arc},samples(end,:));
         L(i,arc) = energies(end); %Note the pseudotimes not updates so prior prob. the same.
         Param.Store{count} = Data; %Store every 100th sample
         clear global HMC_MOM
         end
         disp(['Step ' num2str(i)])
         end        
         
         %Save every 10 steps
         if double(int64(i/10))==(i/10)
             Data.update.uc1 = uc1;
             Data.update.uc2 = uc2;
             Data.update.uc3 = uc3;
             Data.update.uc4 = uc4;
             Data.update.uc5 = uc5;
             Data.update.uc6 = uc6;
             Data.update.uc7 = uc7;
             Data.update.ub1 = ub1;
             Data.update.ub4 = ub4;
             Data.update.count = count;
             Output.Data  = Data;
             Output.L     = L;
             Output.Pr    = Prt;
             Output.Param = Param;
             Output.Xstar = Xstar;
             Output.stepno = i;
             Output.s = rng;
            
            save(['Marker_Pseudotime_' num2str(seed) '_' num2str(missingT) '_cov2E_pr.mat'],'Output')            
            %save(['v4Marker_Pseudotime_percent_run=' num2str(seed) '_withESC_TimeGaussUpdate_TimeSwapUpdate_BranchUpdate_withHyperparmsII_withoutBulk_missingt6_cov5a_extraupdates_transformdropout.mat'],'Output')
         end
         
         
end


            

Output.Data  = Data;
Output.L     = L;
Output.Pr    = Prt;
Output.Param = Param;
Output.Xstar = Xstar;
Output.stepno = i;
Output.s = rng;

save(['Marker_Pseudotime_' num2str(seed) '_' num2str(missingT) '_cov2E_pr.mat'],'Output')            
%save(['Marker_Pseudotime_' num2str(seed) '_' num2str(missingT) '_cov5a.mat'],'Output')            
%save(['v4Marker_Pseudotime_percent_run=' num2str(seed) '_withESC_TimeGaussUpdate_TimeSwapUpdate_BranchUpdate_withHyperparmsII_withoutBulk_missingt6_cov5c_extraupdates_transformdropout.mat'],'Output')


function [Data,L,Prt,updated] = updateTimeSwap(Data,Lin,Prtin);

FL      = Data.update.FixTime;
Times   = Data.update.t;
uT      = unique(Times);

X = Data.update.x;
Y = Data.update.y;

Set1 = find(X(:,2)==1 & FL==0);
Set2 = find(X(:,2)==2 & FL==0);

%Make the update one branch or the other, not both.

if length(Set1)>=2 %Update time of branch 1    
    ind = randperm(length(Set1));
    if Set1(ind(1))==min(Set1)
        X([Set1(1),Set1(2)],1) = X([Set1(2),Set1(1)],1);
        Rat1 = 1;
    elseif Set1(ind(1))==max(Set1)
        X([Set1(end),Set1(end-1)],1) = X([Set1(end-1),Set1(end)],1);   
        Rat1 = 1;
    else
        choice = sign(randi([0,1])-0.5);
        indalt=find(Set1==Set1(ind(1)));        
        X([Set1(ind(1)),Set1(indalt+choice)],1) = X([Set1(indalt+choice),Set1(ind(1))],1);
        Rat1 = 0.5;
    end
else
    Rat1 = 1;
end


if length(Set2)>=2 %Update time for branch 2
    ind = randperm(length(Set2));
    if Set2(ind(1))==min(Set2)
        X([Set2(1),Set2(2)],1) = X([Set2(2),Set2(1)],1);
        Rat2 = 1;
    elseif Set2(ind(1))==max(Set2)
        X([Set2(end),Set2(end-1)],1) = X([Set2(end-1),Set2(end)],1);        
        Rat2 = 1;
    else
        choice = sign(randi([0,1])-0.5);
        indalt=find(Set2==Set2(ind(1)));        
        X([Set2(ind(1)),Set2(indalt+choice)],1) = X([Set2(indalt+choice),Set2(ind(1))],1);
       Rat2 = 0.5;
    end  
else
    Rat2 = 1;
end
    L2 = zeros(1,size(Y,2));
    for j = 1:size(Y,2)        
        ind = find(Y(:,j)~=0);
        par = {'meanConst','covBranchingProcess_2E','likGauss', X(ind,:),Y(ind,j)};
        [L2(1,j) dL] = gp(Data.update.hyp{j}, Data.update.im, par{:}); 
    end
    
    Prt = Prtin;
    for i= 1:length(Data.update.origt)
        if  Data.update.origt(i)==min(Data.update.origt)%Data.update.origt(index1(i))==min(Data.update.origt) %Unlabelled time point
            Prt(1,i) = feval(Data.update.p1{:},X(i,1));
        else
            pri = {@priorGauss,Data.update.origt(i,1),.01};
            Prt(1,i) = feval(pri{:},X(i,1));
        end    
    end    
        
    AR = min(log(1),-sum(L2)+sum(Prt)+sum(Lin)-sum(Prtin)+log(Rat1)+log(Rat2));
    if log(rand(1,1))<AR %Accept
        Data.update.x = X;
        updated = 1;        
        L = L2;
    else
        L = Lin;
        updated = 0;
        Prt = Prtin;
    end
return

function [Data,L,updated] = updateGlobalAssignment(Data,Lin);

FL      = Data.update.FixLabel;
Times   = Data.update.t;
uT      = unique(Times);

X = Data.update.x;
Y = Data.update.y;

if rand(1,1)<0.5
Set1 = find(X(:,2)==1);
inds1 = Set1(randperm(length(Set1)));
N    = binornd(length(inds1)-1,Data.update.pgs)+1; %Number to update
X(inds1(1:N)) = 2;
FT = binopdf(N,length(inds1)-1,Data.update.pgs);
BT = binopdf(N,length(inds1)-1+N,Data.update.pgs);
else
Set1 = find(X(:,2)==2);
inds1 = Set1(randperm(length(Set1)));
N    = binornd(length(inds1)-1,Data.update.pgs)+1; %Number to update 
X(inds1(1:N)) = 1;
FT = binopdf(N,length(inds1)-1,Data.update.pgs);
BT = binopdf(N,length(inds1)-1+N,Data.update.pgs);
end

    L2 = zeros(1,size(Y,2));
    for j = 1:size(Y,2)                
        ind = find(Y(:,j)~=0);
        par = {'meanConst','covBranchingProcess_2E','likGauss', X(ind,:),Y(ind,j)};
        [L2(1,j) dL] = gp(Data.update.hyp{j}, Data.update.im, par{:}); 
    end
            
    AR = min(log(1),-sum(L2)+sum(Lin)+log(BT)-log(FT));
    if log(rand(1,1))<AR %Accept
        Data.update.x = X;
        updated = 1;        
        L = L2;
    else
        L = Lin;
        updated = 0;
    end    
return


function [Data,L,updated] = updateOutliers(Data,Lin,Prtin);

FL      = Data.update.FixLabel;
Times   = Data.update.t;
uT      = unique(Times);

X = Data.update.x;
Y = Data.update.y;

    D = [];
    for j = 1:size(Y,2)     
        
        ind = find(Y(:,j)~=0);
        par = {'meanConst','covBranchingProcess_2E','likGauss', X(ind,:),Y(ind,j),X(:,:)};
        [ymu ys2 fmu fs2] = gp(Data.update.hyp{j}, Data.update.im, par{:});         
        D = [D,fmu];
    end
    Dist = 1./sum(D,2);
    Dist = Dist./sum(Dist);
    
    %Try updating outliers
    ind = gendist(Dist',1,1);
    
    %try flipping this (can generalise and flip multuple outliers).
    if X(ind,2)==1
    X(ind,2) = 1;
    else
    X(ind,2) = 2;
    end
    
    FT = Dist(ind);
    D2 = [];
    for j = 1:size(Y,2)        
        
        ind = find(Y(:,j)~=0);
        par = {'meanConst','covBranchingProcess_2E','likGauss', X(ind,:),Y(ind,j),X(:,:)};
        [ymu ys2 fmu fs2] = gp(Data.update.hyp{j}, Data.update.im, par{:});         
        D2 = [D2,fmu];
    end
    Dist2 = 1./sum(D2,2);
    Dist2 = Dist2./sum(Dist2);
    
    BT = Dist2(ind);
       
    L2 = zeros(1,size(Y,2));
    for j = 1:size(Y,2)        
        
        ind = find(Y(:,j)~=0);
        par = {'meanConst','covBranchingProcess_2E','likGauss', X(ind,:),Y(ind,j)};
        [L2(1,j) dL] = gp(Data.update.hyp{j}, Data.update.im, par{:}); 
    end    
    
    for i= 1:length(Data.update.origt)
        if  Data.update.origt(i)==min(Data.update.origt)%Unlabelled time point
            Prt(1,i) = feval(Data.update.p1{:},X(i,1));
        else
            pri = {@priorGauss,Data.update.origt(i,1),.01};
            Prt(1,i) = feval(pri{:},X(i,1));
        end    
    end    
        
    AR = min(log(1),-sum(L2)+sum(Lin)+log(BT)-log(FT));
    if log(rand(1,1))<AR %Accept
        Data.update.x = X;
        updated = 1;        
        L = L2;
    else
        L = Lin;
        updated = 0;
        Prt = Prtin;
    end   

return

function [Data,L,Prt,updated] = updateGlobalTimeSwap(Data,Lin,Prtin);

FL      = Data.update.FixTime;
Times   = Data.update.t;
uT      = unique(Times);

X = Data.update.x;
Y = Data.update.y;

Set1 = find(FL==0);

%Make the update one branch or the other, not both.
    inds1 = (randperm(length(Set1)));
    X([Set1(inds1(1)),Set1(inds1(2))],1) = X([Set1(inds1(2)),Set1(inds1(1))],1);

    L2 = zeros(1,size(Y,2));
    for j = 1:size(Y,2)        
        
        ind = find(Y(:,j)~=0);
        par = {'meanConst','covBranchingProcess_2E','likGauss', X(ind,:),Y(ind,j)};
        [L2(1,j) dL] = gp(Data.update.hyp{j}, Data.update.im, par{:}); 
    end
    
    Prt = Prtin;
    for i= [Set1(inds1(1)),Set1(inds1(2))]%1:length(Data.update.origt)
        if  Data.update.origt(i)==min(Data.update.origt)%Data.update.origt(index1(i))==min(Data.update.origt) %Unlabelled time point
            Prt(1,i) = feval(Data.update.p1{:},X(i,1));
        else
            pri = {@priorGauss,Data.update.origt(i,1),.01};
            Prt(1,i) = feval(pri{:},X(i,1));
        end    
    end    
        
    AR = min(log(1),-sum(L2)+sum(Prt)+sum(Lin)-sum(Prtin));
    if log(rand(1,1))<AR %Accept
        Data.update.x = X;
        updated = 1;        
        L = L2;
    else
        L = Lin;
        updated = 0;
        Prt = Prtin;
    end
    

return

function [Data,L,Prt,updated] = updateGlobalTimeGauss(Data,Lin,Prtin);

FL  = Data.update.FixTime; %Which cells are okay to update?
ind = find(FL==0);

%Random subset of cells
ind2 = ind(randperm(length(ind))); %Select these at random to update

%Apply Gaussian perturbation of the subset
X = Data.update.x;
X(ind2(1),1) = X(ind2(1),1) + randn(1,1)*Data.update.TR2;
Y = Data.update.y; 

%Get likelihoods
L2 = zeros(1,size(Y,2));
for arc = 1:size(Y,2)        
    
        ind = find(Y(:,arc)~=0);
        par = {'meanConst','covBranchingProcess_2E','likGauss', X(ind,:),Y(ind,arc)};
        [L2(1,arc) dL] = gp(Data.update.hyp{arc}, Data.update.im, par{:}); 
end
    
Prt = Prtin;
%Get pseudotime prior probs. Vectorise?
for arc = ind2(1)
        if  Data.update.origt(arc)==min(Data.update.origt)%ESCs have sepearate prior
            Prt(1,arc) = feval(Data.update.p1{:},X(arc,1));
        else
            pri = {@priorGauss,Data.update.origt(arc,1),.01}; %All other cells are Gaussian dist.
            Prt(1,arc) = feval(pri{:},X(arc,1));
        end    
end   
        
    deltaE = min(log(1),sum(-L2)+sum(Lin)+sum(Prt)-sum(Prtin));
    if log(rand(1,1))<deltaE %Accept (update the ordering)      
        Data.update.x(ind2(1),1) = X(ind2(1),1);            
        updated = 1;        
        L = L2;
    else
        %Reject the update (can optimise speed later)
        L = Lin;
        updated = 0;
        Prt = Prtin;
    end
    
return


function [Data,L,Prt,updated] = updateTimeGauss(Data,Lin,Prtin);

FL  = Data.update.FixTime; %Which cells are okay to update?
ind = find(FL==0);

%Random subset of cells
N    = binornd(length(ind)-1,Data.update.p)+1; %Number to update
ind2 = ind(randperm(length(ind))); %Select these at random to update

%Apply Gaussian perturbation of the subset
X = Data.update.x;
X(ind2(1:N),1) = X(ind2(1:N),1) + randn(N,1)*Data.update.TR;
Y = Data.update.y; 

%Get likelihoods
L2 = zeros(1,size(Y,2));
for arc = 1:size(Y,2)        
    
        ind = find(Y(:,arc)~=0);
        par = {'meanConst','covBranchingProcess_2E','likGauss', X(ind,:),Y(ind,arc)};
        [L2(1,arc) dL] = gp(Data.update.hyp{arc}, Data.update.im, par{:}); 
end  
    
Prt = Prtin;
%Get pseudotime prior probs. Vectorise?
for arc = ind2(1:N)' %1:length(Data.update.origt)
        if  Data.update.origt(arc)==min(Data.update.origt)%Data.update.origt(index1(arc))==min(Data.update.origt) %ESCs have sepearate prior
            Prt(1,arc) = feval(Data.update.p1{:},X(arc,1));
        else
            pri = {@priorGauss,Data.update.origt(arc,1),.01}; %All other cells are Gaussian dist.
            Prt(1,arc) = feval(pri{:},X(arc,1));
        end    
end   
        
    deltaE = min(log(1),sum(-L2)+sum(Lin)+sum(Prt)-sum(Prtin));
    if log(rand(1,1))<deltaE %Accept (update the ordering)      
        Data.update.x(ind2(1:N),1) = X(ind2(1:N),1);          
        updated = 1;        
        L = L2;
    else
        %Reject the update (can optimise speed later)
        L = Lin;
        updated = 0;
        Prt = Prtin;
    end
    
return

function [Data,L,updated] = updateAssignment(Data,Lin);

FL      = Data.update.FixLabel;
Times   = Data.update.t;
uT      = unique(Times);
inds2up = find(FL==0);
u       = inds2up(randi([1 length(inds2up)],1,1)); %Randomly switch assignment of one data point

if Data.update.x(u,2)==1 %Current update is branch 1
        
        %L1 = Lin;            %Current -ve log likelihood for B1
        %X1 = Data.update.x;  %X and Y remain in the same order
        %Y1 = Data.update.y; 

        X2 = Data.update.x; %Try switching branch assignment
        Y2 = Data.update.y; %Update X and Y
        X2(u,2)=2;          %Can we make this bit faster? Only swapping 1 cell.        
    
    for j = 1:size(Y2,2)        
        %ind = find((Y2(:,j))~=pi);
        
        ind = find(Y2(:,j)~=0);
        par = {'meanConst','covBranchingProcess_2E','likGauss', X2(ind,:),Y2(ind,j)};
        [L2(1,j) dL] = gp(Data.update.hyp{j}, Data.update.im, par{:}); 
    end    
    
    scale  = sum(Lin); %Numerical problems so rescale.
    P      = log(exp(-sum(Lin)+scale)./(exp(-sum(Lin)+scale)+exp(-sum(L2)+scale)));
             %Probability of assigning to branch 1    
    choice = log(rand(1,1)); %Randomly assign ...

    if choice<P %Assign to branch 1 i.e., do not update. Likeihood stays same     
        updated = 0;
        L       = Lin;
    else %Assign to branch 2 ...
        Data.update.x(u,2) = 2;% = X2; %Update datasets        
        L = L2;
        updated = 1;        
    end
            
elseif Data.update.x(u,2)==2 %Current locus is Branch 2

    %L2 = Lin;
    %X2 = Data.update.x;
    %Y2 = Data.update.y;
    

    X1 = Data.update.x;
    Y1 = Data.update.y;    
    X1(u,2)=1; %Update to branch 1
    %[X1,index1] = sortrows(X1,[2,1]);
    %Y1 = Y1(index1,:);        
    
    for j = 1:size(Y1,2)
        %ind = find((Y1(:,j))~=pi);
        
        ind = find(Y1(:,j)~=0);
        par = {'meanConst','covBranchingProcess_2E','likGauss', X1(ind,:),Y1(ind,j)};
        [L1(1,j) dL] = gp(Data.update.hyp{j}, Data.update.im, par{:});
    end

    scale  = sum(L1);
    P      = log(exp(-sum(L1)+scale)./(exp(-sum(L1)+scale)+exp(-sum(Lin)+scale)));
    choice = log(rand(1,1));

    if choice<P %Assign to branch 1
        Data.update.x(u,2)=1;% = X1;      
        updated = -1;        
        L = L1;
    else
        updated = 0;
        L       = Lin;
    end    
end


function [Data,Xin,Yin,Xstar] = initialiseBranch(y1,t1,Type,Sex,uID,FixLabel,FixTime,Assign,prior);

origt = t1;
gridt = linspace(0,1,50);
gridLab = [ones(1,50),2*ones(1,50)];
gridX = [gridt,gridt; ones(1,50),2*ones(1,50)]';

p1 = {@priorSmoothBox2,0,1,30};

Data.orig.hyp = cell(1,1);
for arc = 1:size(y1,1)

    %Initialise hyperparams
      %Initialise hyperparams
    hyp.mean = mean(y1(arc,:));
    hyp.lik  = 1;
    hyp.cov  = [ 0.45; 4; 0.5; 1; 0.5; 1];

    %Only take those that are fixed and expressed
    Xtrain = [origt(find(FixLabel==1 & y1(arc,:)~=0));Assign(find(FixLabel==1 & y1(arc,:)~=0))]';          
    Ytrain = y1(arc,find(FixLabel==1 & y1(arc,:)~=0))';      
    %Initilise hyperparams (with a cheeky initial optimisation)
    im  = {@infPrior,@infExact,prior};
    
    par = {'meanConst','covBranchingProcess_2E','likGauss', Xtrain,Ytrain};    
    hyp_init = feval(@minimize, hyp, @gp, -10000, im, par{:});    
    Data.orig.hyp{arc} = hyp_init;    
end

%Now we are going to assign the missing data. Initially this is over a grid
%of timepoints.
for cells = 1:length(FixLabel);%find(FixLabel==0) %These have no times or labels ...  
    
    if FixLabel(cells)==0
    
    Lik = [];

    for arc = 1:size(y1,1)    
    %Put training data in correct format
    Xtrain = [origt(find(FixLabel==1 & y1(arc,:)~=0));Assign(find(FixLabel==1 & y1(arc,:)~=0))]';          
    Ytrain = y1(arc,find(FixLabel==1 & y1(arc,:)~=0))';      
    
    %Prediction points
    if y1(arc,cells)==0 %Missing data point, uniform distribution
        lp = log(ones(2*length(gridt))./2*length(gridt));
        Lik = [Lik;lp'];
    else
        %Xtest = [unique(t1(find(t1>=0)))',ones(length(unique(t1(find(t1>=0)))),1);unique(t1(find(t1>=0)))',2*ones(length(unique(t1(find(t1>=0)))),1)];
        %Ytest = repmat(y1(arc,cells),size(Xtest,1),1); 
        Ytest = repmat(y1(arc,cells),size(gridX,1),1);   

        im  = {@infPrior,@infExact,prior};
        par = {'meanConst','covBranchingProcess_2E','likGauss', Xtrain, Ytrain, gridX, Ytest};    
        [ymu ys2 fmu fs2 lp] = gp(Data.orig.hyp{arc}, im, par{:});     

        Lik = [Lik.^0.1;lp'];    
    end
    end
    
        Prt = zeros(1,2*length(gridt));
        if  origt(cells)==min(origt)%Unlabelled time point
            Prt(1,:) = feval(p1{:},[gridt,gridt]);
        else
            pri = {@priorGauss,origt(cells),.01};
            Prt(1,:) = feval(pri{:},[gridt,gridt]);
        end  
        Lik = [Lik;Prt];
        
    %Now assign time and label to the cell of interest    
    Tim = [gridt,gridt]; %[unique(t1(find(t1>=0))),unique(t1(find(t1>=0)))];
    %gridLab = %[ones(1,length(unique(t1(find(t1>=0))))),2*ones(1,length(unique(t1(find(t1>=0)))))];
    P = exp(sum(Lik,1)-max(sum(Lik,1)))./sum(exp(sum(Lik,1)-max(sum(Lik,1))));%exp(sum(Lik,1))./sum(exp(sum(Lik,1)));     
    ind = gendist(P,1,1);
    %ML or sample from distribution?    
    t1(cells) = Tim(ind);
    Assign(cells) = gridLab(ind);    
    
    else
        
        
 Lik = [];
    
    for arc = 1:size(y1,1)    
    %Put training data in correct format
    Xtrain = [origt(find(FixLabel==1 & y1(arc,:)~=0));Assign(find(FixLabel==1 & y1(arc,:)~=0))]';          
    Ytrain = y1(arc,find(FixLabel==1 & y1(arc,:)~=0))';      
    
    %Prediction points
    if y1(arc,cells)==0 %Missing data point, uniform distribution
        lp = log(ones(length(gridt))./length(gridt));
        Lik = [Lik;lp'];
    else
        %Xtest = [unique(t1(find(t1>=0)))',ones(length(unique(t1(find(t1>=0)))),1);unique(t1(find(t1>=0)))',2*ones(length(unique(t1(find(t1>=0)))),1)];
        Ytest = repmat(y1(arc,cells),size(gridX,1)/2,1);    

        im  = {@infPrior,@infExact,prior};
        
        if Assign(cells)==1
        Xtest=gridX(1:50,:);
        else
        Xtest=gridX(51:100,:);
        end
        par = {'meanConst','covBranchingProcess_2E','likGauss', Xtrain, Ytrain, Xtest, Ytest};    
        [ymu ys2 fmu fs2 lp] = gp(Data.orig.hyp{arc}, im, par{:});     

        Lik = [Lik;lp'];    
        
        %for i= 1:length(Data.update.origt)
        %if  Data.update.origt(arc)==min(Data.update.origt)%Unlabelled time point
        %    Prt(1,i) = feval(Data.update.p1{:},X(i,1));
        %else
        %    pri = {@priorGauss,Data.update.origt(i,1),.01};
        %    Prt(1,i) = feval(pri{:},X(i,1));
        %end    
        %end   
        
        
    end
    end
    
    try

        Prt = zeros(1,length(gridt));
        if  origt(cells)==min(origt)%Unlabelled time point
            Prt(1,:) = feval(p1{:},[gridt]);
        else
            pri = {@priorGauss,origt(cells),.01};
            Prt(1,:) = feval(pri{:},[gridt]);
        end  
        Lik = [Lik.^0.1;Prt]; %Initialise over power posterior
    catch
        keyboard
    end
        
    %Now assign time and label to the cell of interest    
    Tim = [gridt]; %[unique(t1(find(t1>=0))),unique(t1(find(t1>=0)))];
    %gridLab = %[ones(1,length(unique(t1(find(t1>=0))))),2*ones(1,length(unique(t1(find(t1>=0)))))];
    P = exp(sum(Lik,1)-max(sum(Lik,1)))./sum(exp(sum(Lik,1)-max(sum(Lik,1))));%exp(sum(Lik,1))./sum(exp(sum(Lik,1)));     
    ind = gendist(P,1,1);
    %ML or sample from distribution?    
    t1(cells) = Tim(ind);
    %Assign(cells) = gridLab(ind);                         
    end
    
    
    
    
end

%Get the starting likelihood

Xin = [t1',Assign'];
Yin = y1';

Data.orig.t = t1';%(index)';
Data.orig.origt = origt';%(index)';
Data.orig.y = Yin;
Data.orig.Sex = Sex';%(index)';
Data.orig.Type = Type';%(index)';
Data.orig.FixLabel = FixLabel';%(index)';
Data.orig.FixTime = FixTime';%(index)';
Data.orig.Assign = Assign;%(index);
Data.orig.uID = uID';%(index)';%';
Data.orig.x = Xin;
Data.orig.tgrid = unique(t1);

%Information about priors
Data.orig.cov = 'covBranchingProcess_2E';
Data.orig.im = {@infPrior,@infExact,prior};
Data.orig.p1 = p1;%{@priorSmoothBox2,0,1,30};
Data.orig.burnin = 1000; %Burnin
Data.orig.TR = 0.01;  %Standard deviation of steps in time perturbation
Data.orig.TR2 = 0.1;  %Standard deviation of steps in time perturbation
Data.orig.p  = 0.02;  %Probabilty for binomial sampling of cells to update
Data.orig.pgs = 0.02;

Xstar    = [linspace(0,1,100)',ones(1,100)';linspace(0,1,100)',2*ones(1,100)'];
%hyp.mean = mean(X(:,1));
%hyp.lik  = -0.2097;
%hyp.cov  = [ 0.4561; 4.3786; 4.4966; 1.4560; -0.4405; -0.7716];
%Data.orig.hyp = hyp;
Data.update   = Data.orig;

function [options, prior, altprior] = initialiseprior(D1);

%HMC options
options     = foptions; % Default options vector.
options(1)  = 0;		% Switch on diagnostics.
options(5)  = 0;		% Use persistence
options(7)  = 3;	    % Number of steps in trajectory.
options(14) = 1;	    % Number of Monte Carlo samples returned. 
options(15) = 1000;	    % Number of samples omitted at start of chain.
options(18) = 0.02;

%Mean and variance of the means (for prior distributions). Empirically set.
Param.HYPstore = [];
m_mu    = mean(mean(D1.data(4:end,:),2));
s_mu    = var(mean(D1.data(4:end,:),2));
m_noise = 2;      %Noise hyperparamst
s_noise = 1;

m_time  = 0.4; %Switch time hyperparms, centred around 1st stage of PGC (0.5385)
s_time  = 1;

m_switch= 4;      %Switch rate
s_switch= 1;

m_LS1   = 0.5;      %LS 1
s_LS1   = .05;
m_V1    = 1;      %Process variance 1
s_V1    = 1;

%Prior over hyperparameters
p1 = {@priorGauss,m_mu,s_mu};         %Mean
p2 = {@priorGauss,m_noise,s_noise};   %Noise
p3 = {@priorSmoothBox1,m_time,s_time,10};%Switch time, smooth between blastocyst and late. Lower probability of being earlier. 
p4 = {@priorGauss,m_switch,s_switch}; %Switch rate between processes
p5 = {@priorGauss,m_LS1,s_LS1};       %Length-scale 1
p6 = {@priorGauss,m_V1,s_V1};         %Process-variance 1
prior.mean = {p1};                    %Mean hyp.
prior.lik =  {p2};                    %Likelihood hyp.
prior.cov  = {p3;p4;p5;p6;p5;p6};     %Covariance hyp.


p1_alt = {@priorClamped};             %Mean
altprior.mean = {p1_alt};
altprior.lik  = {p1_alt};
altprior.cov  = {p3,p4,p1_alt,p1_alt,p1_alt,p1_alt};

function [D1,t1,Type,Sex,uID,FixLabel,FixTime,Assign,nocells] = initialisedata(missingT);


%Load the data (in this case PGC markers)
%D1    = importdata('/Users/christopher_penfold/Desktop/Lab Meeting/NFKB_IL_WNT.xlsx');
 D1    = importdata('./MarkerGenes_Naokos_Paper.xlsx');
% Data is in the following format:
% (No. genes + 3) x (No. cells) 
%
% Row 1 (Type): ESC/unlabelled (-1) Preimplantation (0) PGC (1) or Soma (2)
% Row 2 (Sex):  Unkown (0) Female (1) Male (2)
% Row 3 (Capture time): ESC (-1), Preimplantation 0-6, specfied cell PGC or soma (>6)
% Row 4 onwards: Gene expression levels

ind_missing = find(D1.data(3,:)==-1); %Get indices for the ESCs/unlabelled. These will be assigned on the fly ...
ind_obs     = setdiff(1:1:size(D1.data,2),ind_missing); %Index for everything else

t1   = D1.data(3,:);            %Vector of developmental stage (capture time)
D1.data(4:end,:) = log2(D1.data(4:end,:)+1); %log_2 transform expression data
%y1   = D1.data(4:end,:);       %Expression levels. Duplicate (remove)
Type = D1.data(1,:);            %Cell type
Sex  = D1.data(2,:);            %Sex
uID  = 1:1:length(t1);          %Unique ID (for rejumbling everything later on)

%inds = find(Sex==0 || Sex==2);

uT = unique(t1);
uType = unique(Type);
uSex = unique(Sex);


% for ii = 1:length(unique(t1))
%     for jj = 1:length(unique(Type))
%         for kk = 1:length(unique(Sex))
%             ind = find(t1==uT(ii) & Sex==uSex(kk) & Type==uType(jj));            
%             for ll = 1:size(D1.data(4:end,:),1)                
%                 %Get rid of extreme dropout
%                 y = D1.data(ll+3,ind);                
%                 if length(find(y~=0))>3 %If there are non zeros                
%                 y(find(y==0)) = median(y(find(y~=0)));
%                 D1.data(ll+3,ind) = y; 
%                 end
%             end
%             
%         end
%     end
% end
     


FixLabel = zeros(1,length(t1)); %Fix labels for a subset of cells where FixLabel=1
FixLabel(1,ind_obs) = ones;     %i.e., do not update the label in the Gibbs step
FixTime = zeros(1,length(t1));  %Fix these time points (will be used for bulk data)
                                %i.e. do not update the time point

%Need to assign branch label (1 or 2). Do this randomly for the pre-implantation.
Assign  = Type;
Assign(find(Assign==0 | Assign==-1)) = randi([1 2],1,length(find(Assign==0 | Assign==-1))); %Random assignment of pre-implantation and ESCs

%Okay, lets pretend we have bulk data (actually this is the average of the single cell data). Maybe a little naughty.
%  Xin = [];
%  for i = 1:14 %There are 14 capture times
%      Dat = D1.data(4:end,:);    
%      if i <= 6
%              ind1 = find(t1==(i-1)); %Get all preimplantation at timepoint (i-1)
%              X    = mean(Dat(:,ind1),2); %Expression level for all genes
%              Xin  = [Xin,[randi([1 2]); 1; i; X]]; %Type, sex and capture time
%      else
%          ind1 = find(t1==(i-1) & Type==1); %PGC
%          ind2 = find(t1==(i-1) & Type==2); %Soma
%          X1 = mean(Dat(:,ind1),2);
%          X2 = mean(Dat(:,ind2),2);            
%          Xin = [Xin,[1;1;i;X1],[2;1;i;X2]];
%      end
%  end
%   Xin = Xin(:,find(isnan(sum(Xin,1))==0));
%Finished generating the bulk measurments

%Type(:) = -1; %Okay now lets develop amnesia over these branch labels (blastocyst).
%t1(:) = -1;
Type(find(t1==6 | t1==missingT))     = -1;
FixLabel(find(t1==6 | t1==missingT)) = 0;
t1(find(t1==6 | t1==missingT))       = -1;

t1       = [t1];
Sex      = [Sex];
Type     = [Type];
nocells  = max(uID);
uID      = [uID]; %Unique labels for each cell
D1.data  = [D1.data];
%y1       = [y1,Xin(4:end,:)];
NoObs    = size(FixLabel,2);

%FixLabel(:) = zeros;
FixLabel    = [FixLabel]; %Fix only the bulk data
FixTime     = [FixTime];  
Assign      = [Assign]; %????? Is this random? Do we use this variable?

%Concatenate bulk data with single cell data
%t1       = [t1,Xin(3,:)];
%Sex      = [Sex,Xin(2,:)];
%Type     = [Type,Xin(1,:)];
%nocells  = max(uID);
%uID      = [uID,nocells+1:size(t1,2)]; %Unique labels for each cell
%D1.data  = [D1.data,Xin];
%%y1       = [y1,Xin(4:end,:)];
%NoObs    = size(FixLabel,2);

%FixLabel(:) = zeros;
%FixLabel    = [FixLabel,ones(1,size(Xin,2))]; %Fix only the bulk data
%FixTime     = [FixTime,ones(1,size(Xin,2))];  
%Assign      = [Assign,Xin(1,:)]; %????? Is this random? Do we use this variable?

t1          = t1./max(t1); %Re-scale time beteween 0 and 1

return
