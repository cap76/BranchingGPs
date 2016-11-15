addpath(genpath('../'))
addpath(genpath('/Users/christopher_penfold/Desktop/Code/gpss-research/source/gpml'))
addpath(genpath('/Users/christopher_penfold/Desktop/Code/deepGP/netlab3_3/'))
warning off all

%Input locations
x = linspace(0,10,100);

%Covariance function (this is the same as covBranchingProcess_2B
covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
covN2 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
K1   = covSEiso([1,3],x');
K2   = feval(covN1{:}, [4,0.5,1.8,1], x');
K3   = feval(covN1{:}, [4,0.5,0.8,1], x');
Kall = [K1,K1;K1,K1+K2];

%Generate noisy training data
X = [x',ones(100,1);x',2*ones(100,1)];
y = real(gsamp(zeros(1,200),Kall,10));
Y = y(1,1:200)' + rand(length(y(1,1:200)),1)*1;

%Prediction
Xstar = [linspace(0,10,1000)',ones(1000,1);linspace(0,10,1000)',2*ones(1000,1)];

%Initialise hyperparameters
hyp.cov = [8;0.5;rand(4,1)];%[3,0.5,0.5,1,3,0.3,2,2,2]';
hyp.mean = mean(X(:,1));
hyp.lik = .1;

%Optimise hyperparameters
hyp_pN      = feval(@minimize, hyp, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_2A','likGauss',X,Y);

%Predictions
[ymu ys2 fmu fs2   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingProcess_2A', 'likGauss', X, Y, Xstar);

%Now plot
K = covBranchingProcess_2A(hyp_pN.cov, X);
subplot(1,3,1); imagesc(K);
y0 = real(gsamp(zeros(1,200),K,1));
subplot(1,3,2); hold on
plot(X(:,1),Y,'ks')
xlim([0 10])
subplot(1,3,3); 
errorbar(Xstar(:,1),ymu,2*sqrt(ys2),'b.')
hold on
xlim([0 10])
title([hyp_pN.cov(1)])