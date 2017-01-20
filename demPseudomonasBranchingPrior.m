function [Fin] = demPseudomonasBranchingPrior(batchi);

%Fit branching process to the Arabidopsis/pseudomonas dataset (Lewis et al., 2015).
%Split the analysis into batches of 1000 for running on the cluster. 

addpath(genpath('../'))
%addpath(genpath('../Pseudotime_Algorithms/Functions/gpml-matlab-v3.6-2015-07-07/'))
%addpath(('./gpml/cov'))
%addpath(genpath('./netlab3_3/'))
 
D1 = importdata('./Pa13-Combo2.txt');

%try %Try to resume analysis 
%    load(['./results/pseudomonas/PseudomonasResults_Oct_5_' num2str(batchi) '.mat'])
%    startind = length(Ps)+1;
%    endind   = (batchi)*1000;
%catch
%    startind = (batchi-1)*1000 + 1;
%    endind   = (batchi)*1000;
%endPseudomonas
%load('/Users/christopher_penfold/Desktop/GitHub/demos/results/pseudomonas/PseudomonasBranchingResults.mat')
%endind = length(Pseudomonas);
%endind = min(endind,size(D1.data,1));

%Measurement times. Permute these for different branching structures. 
tt = [0,2,3,4,6,7,8,10,11,12,14,16,17.5];
X1 = [repmat(tt,1,12); ones(1,52),2*ones(1,52),3*ones(1,52)]';
X2 = [repmat(tt,1,12); ones(1,52),ones(1,52),2*ones(1,52)]';
X3 = [repmat(tt,1,12); ones(1,52),2*ones(1,52),2*ones(1,52)]';

%Prediction time points.
Xstar1 = [repmat(linspace(0,17.5,50),1,3);ones(1,50),2*ones(1,50),3*ones(1,50)]';
Xstar2 = [repmat(linspace(0,17.5,50),1,3);ones(1,50),1*ones(1,50),2*ones(1,50)]';
Xstar3 = [repmat(linspace(0,17.5,50),1,3);ones(1,50),2*ones(1,50),2*ones(1,50)]';
keyboard
for i = 1:length(Pseudomonas)%startind:endind

     if isempty(Pseudomonas{i})
    
         
         keyboard
pcp1     = {@priorGamma,2,2};    %Mean 4, std 8 
pcp1p2   = {@priorGamma,4,2};    %Mean 8, std 16 
pcp2     = {@priorGauss,1,0.5};  %Quick transitions        
pctheta1 = {@priorGauss,3,1};    %Length-scale              
pctheta2 = {@priorGauss,1,1};    %Variance             
pcmean   = {@priorGauss,0,1};    %Zero, unit-variance
pclik    = {@priorGauss,0.1,0.5};
    
Y1 = D1.data(i,:)'; %Mock hrp DC
Y2 = D1.data(i,[1:52,2*52+1:3*52,52+1:2*52])'; %Mock DC hrpA

%Initialise some parameters
l1 = log(3); l2 = log(3); lg = log(3); v1 = log(3); v2 = log(3); vg = log(std(Y1));

%3 branch structure. Mock->hrp->DC
hyp.cov  = [4;1;8;1;l1;v1;l1;v1;l1;v1]; hyp.mean = mean(Y1(:,1)); hyp.lik  = 2;
prior.mean = {[]};  prior.cov  = {pcp1;[];pcp1p2;[];[];[];[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y1};
par1b = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y1,Xstar1};
hyp_pN1 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
[L1 dL1] = feval(@gp,hyp_pN1, im, par1a{:});         % optimise
[ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, im, par1b{:});

%3-branch, with one main two indendetly emerged branches, Mock->hrp, Mock->DC
par1a = {'meanConst','covBranchingProcess_3B','likGauss',X1,Y1};
par1b = {'meanConst','covBranchingProcess_3B','likGauss',X1,Y1,Xstar1};
hyp_pN2 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
[L2 dL2] = feval(@gp,hyp_pN2, im, par1a{:});         % optimise
[ymu2 ys22 fmu2 fs22   ]= feval(@gp,hyp_pN2, im, par1b{:});
%[L2 dL2   ] = gp(hyp_pN2, 'infExact', 'meanConst', 'covBranchingProcess_v4', 'likGauss', X1, Y1);
%[ymu2 ys22 fmu2 fs22   ] = gp(hyp_pN2, 'infExact', 'meanConst', 'covBranchingProcess_v4', 'likGauss', X1, Y1, Xstar1);

%Now do joint modelling of subsets
%1,2->3, Mock,Hrp -> DC
hyp.cov = [4;1;l1;v1;l1;v1];hyp.mean = mean(Y1(:,1));hyp.lik = 2;
prior.mean = {[]};  prior.cov  = {pcp1;[];[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y1};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y1,Xstar2};
hyp_pN3 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
[L3 dL3] = feval(@gp,hyp_pN3, im, par1a{:});         % optimise
[ymu3 ys23 fmu3 fs23   ]= feval(@gp,hyp_pN3, im, par1b{:});
%hyp.cov = [2;0.5;l1;v1;lg;vg];hyp.mean = mean(Y1(:,1));hyp.lik = 2;
%hyp_pN3      = feval(@minimize, hyp, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_v5','likGauss',X2,Y1);
%[L3 dL3   ] = gp(hyp_pN3, 'infExact', 'meanConst', 'covBranchingProcess_v5', 'likGauss', X2, Y1);
%[ymu3 ys23 fmu3 fs23   ] = gp(hyp_pN3, 'infExact', 'meanConst', 'covBranchingProcess_v5', 'likGauss', X2, Y1, Xstar2);
 
%Shift data around
hyp.cov  = [4;1;8;1.5;l1;v1;l1;v1;l1;v1]; hyp.mean = mean(Y1(:,1)); hyp.lik  = 2;
prior.mean = {[]};  prior.cov  = {pcp1;[];pcp1p2;[];[];[];[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y2};
par1b = {'meanConst','covBranchingProcess_3A','likGauss',X1,Y2,Xstar1};
hyp_pN4 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
%hyp_pN3      = feval(@minimize, hyp, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_v4','likGauss',X1,Y1);
[L4 dL4] = feval(@gp,hyp_pN4, im, par1a{:});         % optimise
[ymu4 ys24 fmu4 fs24   ]= feval(@gp,hyp_pN4, im, par1b{:});

%1->2->3,  Mock ->DC->hrp
%hyp_pN4      = feval(@minimize, hyp, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_v3','likGauss',X1,Y2);
%[L4 dL4   ] = gp(hyp_pN4, 'infExact', 'meanConst', 'covBranchingProcess_v3', 'likGauss', X1, Y2);
%[ymu4 ys24 fmu4 fs24   ] = gp(hyp_pN4, 'infExact', 'meanConst', 'covBranchingProcess_v3', 'likGauss', X1, Y2, Xstar1);

%1->2,1->3,  Mock->DC, Mock->hrp (same as 2)
par1a = {'meanConst','covBranchingProcess_3B','likGauss',X1,Y2};
par1b = {'meanConst','covBranchingProcess_3B','likGauss',X1,Y2,Xstar1};
hyp_pN5 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
%hyp_pN3      = feval(@minimize, hyp, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_v4','likGauss',X1,Y1);
[L5 dL5] = feval(@gp,hyp_pN5, im, par1a{:});         % optimise
[ymu5 ys25 fmu5 fs25   ]= feval(@gp,hyp_pN5, im, par1b{:});
%hyp_pN5      = feval(@minimize, hyp, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_v4','likGauss',X1,Y2);
%[L5 dL5] = gp(hyp_pN5, 'infExact', 'meanConst', 'covBranchingProcess_v4', 'likGauss', X1, Y2);
%[ymu5 ys25 fmu5 fs25   ] = gp(hyp_pN5, 'infExact', 'meanConst', 'covBranchingProcess_v4', 'likGauss', X1, Y2, Xstar1);


%keyboard
%Now do joint modelling of subsets
%1,2->3,   Mock,hrp -> DC
hyp.cov = [4;1;l1;v1;l1;v1];hyp.mean = mean(Y1(:,1));hyp.lik = 2;
prior.mean = {[]};  prior.cov  = {pcp1;[];[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y2,Xstar2};
hyp_pN6 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
[L6 dL6] = feval(@gp,hyp_pN6, im, par1a{:});         % optimise
[ymu6 ys26 fmu6 fs26   ]= feval(@gp,hyp_pN6, im, par1b{:});
%hyp_pN6      = feval(@minimize, hyp, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_v5','likGauss',X2,Y2);
%[L6 dL6] = gp(hyp_pN6, 'infExact', 'meanConst', 'covBranchingProcess_v5', 'likGauss', X2, Y2);
%[ymu6 ys26 fmu6 fs26   ] = gp(hyp_pN6, 'infExact', 'meanConst', 'covBranchingProcess_v5', 'likGauss', X2, Y2, Xstar2);

% Mock-> DC,hrp
hyp.cov = [4;1;l1;v1;l1;v1];hyp.mean = mean(Y1(:,1));hyp.lik = 2;
par1a = {'meanConst','covBranchingProcess_2A','likGauss',X3,Y1};
par1b = {'meanConst','covBranchingProcess_2A','likGauss',X3,Y1,Xstar3};
hyp_pN8 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
[L8 dL8] = feval(@gp,hyp_pN8, im, par1a{:});         % optimise
[ymu8 ys28 fmu8 fs28   ]= feval(@gp,hyp_pN8, im, par1b{:});

%hyp_pN8      = feval(@minimize, hyp, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_v5','likGauss',X3,Y1);
%[L8 dL8] = gp(hyp_pN8, 'infExact', 'meanConst', 'covBranchingProcess_v5', 'likGauss', X3, Y1);
%[ymu8 ys28 fmu8 fs28   ] = gp(hyp_pN6, 'infExact', 'meanConst', 'covBranchingProcess_v5', 'likGauss', X3, Y1, Xstar3);

%Joint GP (not DE). 
hyp.cov = [l1;v1]; hyp.mean = mean(Y1(:,1)); hyp.lik = 2;
prior.mean = {[]};  prior.cov  = {[];[]}; prior.lik = {pclik};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covSEiso','likGauss',X1(:,1),Y1};
par1b = {'meanConst','covSEiso','likGauss',X1(:,1),Y1,Xstar1(:,1)};
hyp_pN7 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
[L7 dL7] = feval(@gp,hyp_pN7, im, par1a{:});         % optimise
[ymu7 ys27 fmu7 fs27   ]= feval(@gp,hyp_pN7, im, par1b{:});

%Store likelihoods
L = -[L1,L2,L3,L4,L5,L6,L7,L8];
AIC = 2*[12,12,8,12,12,8,4,8] - 2*L;
BIC = - 2*L + [12,12,8,12,12,8,4,8]*log(size(X1,1));

Output.L = L;
Output.AIC = AIC;
Output.BIC = BIC;
Output.H1 = hyp_pN1;
Output.H2 = hyp_pN2;
Output.H3 = hyp_pN3;
Output.H4 = hyp_pN4;
Output.H5 = hyp_pN5;
Output.H6 = hyp_pN6;
Output.H7 = hyp_pN7;
Output.H8 = hyp_pN8;
Output.fmu1 = fmu1;
Output.fmu2 = fmu2;
Output.fmu3 = fmu3;
Output.fmu4 = fmu4;
Output.fmu5 = fmu5;
Output.fmu6 = fmu6;
Output.fmu7 = fmu7;
Output.fmu8 = fmu8;
Output.fs21 = fs21;
Output.fs22 = fs22;
Output.fs23 = fs23;
Output.fs24 = fs24;
Output.fs25 = fs25;
Output.fs26 = fs26;
Output.fs27 = fs27;
Output.fs28 = fs28;

Ps{i,1} = Output;
save('/Users/christopher_penfold/Desktop/GitHub/demos/results/pseudomonas/PseudomonasBranchingResults_complete.mat','Pseudomonas')

%Save every 10 steps.

    %save(['./results/pseudomonas/PseudomonasResults_Oct_5_' num2str(batchi) '.mat'],'Ps')

     end

end

save('/Users/christopher_penfold/Desktop/GitHub/demos/results/pseudomonas/PseudomonasBranchingResults_complete.mat','Pseudomonas')
%save(['./results/pseudomonas/PseudomonasResults_Oct_5_' num2str(batchi) '.mat'],'Ps')

Fin = 1;
