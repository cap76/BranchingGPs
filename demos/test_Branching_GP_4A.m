function [Dis, Lout] = test_Branching_GP_4A(seed,plt);

addpath(genpath('../'))
warning off all

rng(seed)

%Input locations
x = linspace(0,10,1000);

%Generate noisy training data (from a constrained model).
X  = [x',ones(1000,1);x',2*ones(1000,1);x',3*ones(1000,1);x',4*ones(1000,1)];
X0 = [x',ones(1000,1);x',ones(1000,1);x',2*ones(1000,1);x',2*ones(1000,1)];
Kall = covBranchingProcess_4A([3,0.5,5,0.5,6,0.5,2,1.8,2,1.8,2,1.8,2,1.8,2,1.8,2,1.8,1,3]', X);    
y = real(gsamp(zeros(1,4000),Kall,100));
Y = y(seed,1:4000)' + rand(length(y(seed,1:4000)),1)*1;

%Subsample the data for training a GP model.
X  = X(1:10:end,:);
X0 = X0(1:10:end,:);
Y  = Y(1:10:end,1);

%Set prediction locations
Xstar = [linspace(0,10,1000)',ones(1000,1);linspace(0,10,1000)',2*ones(1000,1);linspace(0,10,1000)',3*ones(1000,1);linspace(0,10,1000)',4*ones(1000,1)];

%Two-component GP test
Xstar0 = [linspace(0,10,1000)',ones(1000,1);linspace(0,10,1000)',ones(1000,1);linspace(0,10,1000)',2*ones(1000,1);linspace(0,10,1000)',2*ones(1000,1)];

%Initialise hyperparameters for simple 2-branch model
hyp0.cov = [6;0.5;rand(7,1)];
hyp0.mean = mean(X(:,1));
hyp0.lik = .1;
%Optimise hyperparameters
hyp_pN0      = feval(@minimize, hyp0, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_2C','likGauss',X0,Y);
%Predictions
[ymu0 ys20 fmu0 fs20   ] = gp(hyp_pN0, 'infExact', 'meanConst', 'covBranchingProcess_2C', 'likGauss', X0, Y, Xstar0);
[L0 dL  ] = gp(hyp_pN0, 'infExact', 'meanConst', 'covBranchingProcess_2C', 'likGauss', X0, Y);

%Run a joint GP (over union of datasets)
k1={'covMaterniso',3};
hyp1.cov = [rand(2,1)];
hyp1.mean = mean(X(find(X(:,2)==1),1));
hyp1.lik = .1;
%Optimise hyperparameters
hyp_pN1      = feval(@minimize, hyp1, @gp, -10000, @infExact, 'meanConst',k1,'likGauss',X0(:,1),Y);
%Predictions
[ymu1 ys21 fmu1 fs21   ] = gp(hyp_pN1, 'infExact', 'meanConst', k1, 'likGauss', X0(:,1), Y, Xstar0(:,1));
[L1 dL   ] = gp(hyp_pN1, 'infExact', 'meanConst', k1, 'likGauss', X0(:,1), Y);

%Individual GPs of each branch
k1={'covMaterniso',3};
hyp1.cov = [rand(2,1)];
hyp1.mean = mean(X(find(X(:,2)==1),1));
hyp1.lik = .1;
%Optimise hyperparameters
hyp_pN1      = feval(@minimize, hyp1, @gp, -10000, @infExact, 'meanConst',k1,'likGauss',X(find(X(:,2)==1),1),Y(find(X(:,2)==1)));
%Predictions
[ymu1A ys21A fmu1A fs21A   ] = gp(hyp_pN1, 'infExact', 'meanConst', k1, 'likGauss', X(find(X(:,2)==1),1), Y(find(X(:,2)==1)), Xstar(find(Xstar(:,2)==1)));
[L1A dL1   ] = gp(hyp_pN1, 'infExact', 'meanConst', k1, 'likGauss', X(find(X(:,2)==1),1), Y(find(X(:,2)==1)));

hyp1.cov = [rand(2,1)];
hyp1.mean = mean(X(:,1));
hyp1.lik = .1;
%Optimise hyperparameters
hyp_pN1      = feval(@minimize, hyp1, @gp, -10000, @infExact, 'meanConst',k1,'likGauss',X(find(X(:,2)==2),1),Y(find(X(:,2)==2)));
%Predictions
[ymu1B ys21B fmu1B fs21B   ] = gp(hyp_pN1, 'infExact', 'meanConst', k1, 'likGauss', X(find(X(:,2)==2),1), Y(find(X(:,2)==2)), Xstar(find(Xstar(:,2)==2)));
[L1B dL   ] = gp(hyp_pN1, 'infExact', 'meanConst', k1, 'likGauss', X(find(X(:,2)==2),1), Y(find(X(:,2)==2)));
hyp1.cov = [rand(2,1)];
hyp1.mean = mean(X(find(X(:,2)==3),1));
hyp1.lik = .1;
%Optimise hyperparameters
hyp_pN1      = feval(@minimize, hyp1, @gp, -10000, @infExact, 'meanConst',k1,'likGauss',X(find(X(:,2)==3),1),Y(find(X(:,2)==3)));
%Predictions
[ymu1C ys21C fmu1C fs21C   ] = gp(hyp_pN1, 'infExact', 'meanConst', k1, 'likGauss', X(find(X(:,2)==3),1), Y(find(X(:,2)==3)), Xstar(find(Xstar(:,2)==3)));
[L1C dL   ] = gp(hyp_pN1, 'infExact', 'meanConst', k1, 'likGauss', X(find(X(:,2)==3),1), Y(find(X(:,2)==3)));
hyp1.cov = [rand(2,1)];
hyp1.mean = mean(X(find(X(:,2)==4),1));
hyp1.lik = .1;
%Optimise hyperparameters
hyp_pN1      = feval(@minimize, hyp1, @gp, -10000, @infExact, 'meanConst',k1,'likGauss',X(find(X(:,2)==4),1),Y(find(X(:,2)==4)));
%Predictions
[ymu1D ys21D fmu1D fs21D   ] = gp(hyp_pN1, 'infExact', 'meanConst', k1, 'likGauss', X(find(X(:,2)==4),1), Y(find(X(:,2)==4)), Xstar(find(Xstar(:,2)==4)));
[L1D dL   ] = gp(hyp_pN1, 'infExact', 'meanConst', k1, 'likGauss', X(find(X(:,2)==4),1), Y(find(X(:,2)==4)));


%Now do a 4-component GP.
%Initialise hyperparameters
hyp.cov = [6;0.5;6;0.5;7;0.5;7;0.5;6;0.5;6;0.5;rand(14,1)];%[3,0.5,0.5,1,3,0.3,2,2,2]';
hyp.mean = mean(X(:,1));
hyp.lik = .1;

%Optimise hyperparameters
hyp_pN      = feval(@minimize, hyp, @gp, -50000, @infExact, 'meanConst','covBranchingProcess_4B','likGauss',X,Y);

%Predictions
[ymu ys2 fmu fs2   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingProcess_4B', 'likGauss', X, Y, Xstar);
[L dL   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingProcess_4B', 'likGauss', X, Y); 

%Now plot
if plt==1
 K = covBranchingProcess_4B(hyp_pN.cov, X);
 subplot(2,4,1); imagesc(K);
 y0 = real(gsamp(zeros(1,400),K,1));
 subplot(2,4,2); hold on
 plot(X(:,1),Y,'ks')
 xlim([0 10])
 title(['3, 3, 5, 6'])
 subplot(2,4,3); 
 errorbar(Xstar(:,1),ymu,2*sqrt(ys2),'b.')
 hold on
 xlim([0 10])
 title([num2str(hyp_pN.cov(5)) ', ' num2str(hyp_pN.cov(7)) ', ' num2str(hyp_pN.cov(9)) ', ' num2str(hyp_pN.cov(11))])
 
 subplot(2,4,4); 
 errorbar(Xstar0(:,1),ymu0,2*sqrt(ys20),'b.')
 hold on
 xlim([0 10])
 title(['2-component GP'])
  
 subplot(2,4,5); 
 errorbar(Xstar(:,1),ymu1,2*sqrt(ys21),'b.')
 hold on
 xlim([0 10])
 title(['Joint GP'])
  
 subplot(2,4,6); 
 errorbar(Xstar(:,1),[ymu1A;ymu1B;ymu1C;ymu1D],2*sqrt([ys21A;ys21B;ys21C;ys21D]),'b.')
 hold on
 xlim([0 10])
 title(['Independent GP'])
end

%Now get SSE over full range of testpoints
Dis1 = sum((y(seed,:)' - fmu).^2);
Dis2 = sum((y(seed,:)' - fmu0).^2);
Dis3 = sum((y(seed,:)' - fmu1).^2);
Dis4 = sum((y(seed,:)' - [fmu1A;fmu1B;fmu1C;fmu1D]).^2);
Dis = [Dis1,Dis2,Dis3,Dis4];
Lout = [L,L0,L1,L1A,L1B,L1C,L1D];

%Also plot SSE.
if plt==1
subplot(2,4,7); plot([1,2,3,4],[Dis1,Dis2,Dis3,Dis4],'o'), ylabel('SSE')
subplot(2,4,7); plot([1,2,3,4],-[L,L0,L1,L1A+L1B+L1C+L1D],'o'),ylabel('ML')
end