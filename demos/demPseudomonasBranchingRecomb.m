function [Fin] = demPseudomonasBranchingRecomb(batchi);

%addpath(genpath('../Pseudotime_Algorithms/Functions/gpml-matlab-v3.6-2015-07-07/'))
%addpath(('./gpml/cov'))
%addpath(genpath('./netlab3_3/'))
addpath(genpath('../')) 

D1 = importdata('./Pa13-Combo2.txt');

try %Try to resume analysis 
    load(['./results/pseudomonas/PseudomonasResultsRecomb_Feb_5_' num2str(batchi) '_prior1_v2.mat'])
    startind = length(Ps)+1;
    endind   = (batchi)*1000;
catch
    startind = (batchi-1)*1000 + 1;
    endind   = (batchi)*1000;
end
%load('/Users/christopher_penfold/Desktop/GitHub/demos/results/pseudomonas/PseudomonasRecombinationResults.mat')
%load('/Users/christopher_penfold/Desktop/GitHub/demos/results/pseudomonas/PseudomonasBranchingResults.mat')

%endind = length(PseudomonasRecomb);
endind = min(endind,size(D1.data,1));
tt = [0,2,3,4,6,7,8,10,11,12,14,16,17.5];
X1 = [repmat(tt,1,8); ones(1,52),2*ones(1,52)]';

Xstar1 = [repmat(linspace(0,17.5,50),1,2);ones(1,50),2*ones(1,50)]';

for i = startind:endind

    i
    %if isempty(PseudomonasRecomb{i})
    
pcp1     = {@priorGamma,2,2};    % Gaussian prior
pcp1p2   = {@priorGamma,4,2};    %Mean 8, std 16
pcp2_2 = {@priorClamped};        % Gaussian prior
pcp2     = {@priorGauss,0,0.5};  %Quick transitions 

Y2 = D1.data(i,[1:52,2*52+1:3*52,52+1:2*52])'; %Mock, DC, hrpA
Y1 = Y2(1:2*52,:); %Take mock and DC only

l1 = log(3); l2 = log(3); lg = log(3); v1 = log(3); v2 = log(3); vg = log(std(Y1));


%Joint GP
hyp.cov = [lg;vg]; hyp.mean = mean(Y1(:,1)); hyp.lik = 2;
prior.mean = {[]};  prior.cov  = {[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior}; 
par1a = {'meanConst','covSEiso','likGauss',X1(:,1),Y1};
par1b = {'meanConst','covSEiso','likGauss',X1(:,1),Y1,Xstar1(:,1)};
hyp_pN7      = feval(@minimize, hyp, @gp, -40000, im, par1a{:});
[L7 dL7] = gp(hyp_pN7, im, par1a{:});
[ymu7 ys27 fmu7 fs27   ] = gp(hyp_pN7, im, par1b{:});


%Branching process (recombines at 800, essentially at infinity)
hyp.cov  = [800;0; 4;0;hyp_pN7.cov;hyp_pN7.cov]; hyp.mean = mean(Y1(:,1)); hyp.lik  = hyp_pN7.lik; %hyp.cov  = [800;1;1; 4;1;1;l1;v1;l1;v1;l1;v1]; hyp.mean = mean(Y1(:,1)); hyp.lik  = 2;
prior.mean = {[]};  prior.cov  = {pcp2_2;pcp2;pcp1;pcp2;[];[];[];[]}; prior.lik = {[]};%prior.mean = {[]};  prior.cov  = {pcp2_2;pcp2_2;pcp2_2;pcp1;[];[];[];[];[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingRecombinationProcess_2C','likGauss',X1,Y1};
par1b = {'meanConst','covBranchingRecombinationProcess_2C','likGauss',X1,Y1,Xstar1};
hyp_pN3 = feval(@minimize, hyp, @gp, -40000, im, par1a{:});         % optimise
[L3 dL1] = gp(hyp_pN3, im, par1a{:});         % optimise
[ymu3 ys23 fmu3 fs23   ]= feval(@gp,hyp_pN3, im, par1b{:});

%[18;1;3.5;1;l1;v1;l1;v1];
hyp.cov  = [18;0;hyp_pN3.cov(3);0;hyp_pN3.cov(5:end)]; hyp.mean = mean(Y1(:,1)); hyp.lik  = hyp_pN7.lik;
%Branching/recombination process
prior.mean = {[]};  prior.cov  = {pcp1p2;pcp2;pcp1;pcp2;[];[];[];[]}; prior.lik = {[]}; %prior.cov  = {pcp1p2;[];[];pcp1;[];[];[];[];[];[];[];[]}; prior.lik = {[]};
im = {@infPrior,@infExact,prior};                % inference method
par1a = {'meanConst','covBranchingRecombinationProcess_2C','likGauss',X1,Y1};
par1b = {'meanConst','covBranchingRecombinationProcess_2C','likGauss',X1,Y1,Xstar1};
hyp_pN1 = feval(@minimize, hyp, @gp, -40000, im, par1a{:});         % optimise
[L1 dL1] = gp(hyp_pN1, im, par1a{:});         % optimise
[ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, im, par1b{:});



%Get MLs
L = -[L1,L3,L7];
AIC = 2*[10,8,4] - 2*L;
BIC = - 2*L + [10,8,4]*log(size(X1,1));


Output.L = L;
Output.AIC = AIC;
Output.BIC = BIC;
Output.H1 = hyp_pN1;
Output.H3 = hyp_pN3;
Output.H7 = hyp_pN7;
Output.fmu1 = fmu1;
Output.fmu3 = fmu3;
Output.fmu7 = fmu7;
Output.fs21 = fs21;
Output.fs23 = fs23;
Output.fs27 = fs27;

Ps{i,1} = Output;
%Ps{i,1} = Output;

%save(['./results/pseudomonas/PseudomonasResultsRecomb_Oct_5_' num2str(batchi) '.mat'],'Ps')
%save('/Users/christopher_penfold/Desktop/GitHub/demos/results/pseudomonas/PseudomonasRecombinationResults_complete.mat','PseudomonasRecomb')
 save(['./results/pseudomonas/PseudomonasResultsRecomb_Feb_5_' num2str(batchi) '_prior1_v2.mat'],'Ps')

end

 save(['./results/pseudomonas/PseudomonasResultsRecomb_Feb_5_' num2str(batchi) '_prior1_v2.mat'],'Ps')
%save('/Users/christopher_penfold/Desktop/GitHub/demos/results/pseudomonas/PseudomonasRecombinationResults_complete.mat','PseudomonasRecomb')
%save(['./results/pseudomonas/PseudomonasResultsRecomb_Oct_5_' num2str(batchi) '.mat'],'Ps')

Fin = 1;
