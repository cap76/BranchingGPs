addpath(genpath('../'))
addpath(genpath('/Users/christopher_penfold/Desktop/Code/gpss-research/source/gpml'))
addpath(genpath('/Users/christopher_penfold/Desktop/Code/deepGP/netlab3_3/'))
warning off all

%Input locations
x = linspace(0,10,100);

%Covariance function (this is the same as covBranchingProcess_2B
covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
K   = covSEiso([1,3],x');
Kl1   = feval(covN1{:}, [2,0.5,1.8,1], x');
Kl2   = feval(covN1{:}, [2,0.5,1.8,1], x');
K1   = feval(covN1{:}, [4,0.5,1.8,1], x');
K2   = feval(covN1{:}, [4,0.5,0.8,1], x');
K3   = feval(covN1{:}, [5,0.5,0.8,1], x');
K4   = feval(covN1{:}, [5,0.5,0.8,1], x');

Kall = [K+Kl1+K1,K+Kl1,K,K;
        K+Kl1,K+Kl1+K2,K,K;
        K,K,K+Kl2+K3,K+Kl2;
        K,K,K+Kl2,K+Kl2+K4];
    
%Generate noisy training data
X = [x',ones(100,1);x',2*ones(100,1);x',3*ones(100,1);x',4*ones(100,1)];
y = real(gsamp(zeros(1,400),Kall,10));
Y = y(1,1:400)' + rand(length(y(1,1:400)),1)*1;

%Prediction
Xstar = [linspace(0,10,1000)',ones(1000,1);linspace(0,10,1000)',2*ones(1000,1);linspace(0,10,1000)',4*ones(1000,1);linspace(0,10,1000)',4*ones(1000,1)];

%Initialise hyperparameters
hyp.cov = [8;0.5;8;0.5;8;0.5;rand(14,1)];%[3,0.5,0.5,1,3,0.3,2,2,2]';
hyp.mean = mean(X(:,1));
hyp.lik = .1;

%Optimise hyperparameters
hyp_pN      = feval(@minimize, hyp, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_4A','likGauss',X,Y);

%Predictions
[ymu ys2 fmu fs2   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingProcess_4A', 'likGauss', X, Y, Xstar);

%Now plot
K = covBranchingProcess_4A(hyp_pN.cov, X);
subplot(1,3,1); imagesc(K);
y0 = real(gsamp(zeros(1,400),K,1));
subplot(1,3,2); hold on
plot(X(:,1),Y,'ks')
xlim([0 10])
title(['4, 6, 7'])
subplot(1,3,3); 
errorbar(Xstar(:,1),ymu,2*sqrt(ys2),'b.')
hold on
xlim([0 10])
title([num2str(hyp_pN.cov(1)) ', ' num2str(hyp_pN.cov(3)) ', ' num2str(hyp_pN.cov(5))])