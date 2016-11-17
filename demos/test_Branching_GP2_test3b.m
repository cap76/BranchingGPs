function [Dist Lout] = test_Branching_GP2_test3b(seed,plt);

%Asses accuracy of B/RGPs on a correlated dataset with no branching when 
%branch point hyperparameters are fixed at +/- large number i.e., the limit
%B/RGP -> HGP.

rng(seed)

addpath(genpath('../'))

%thetan = log(0.1);

xreal = linspace(-8,8,3000);

N = randperm(3000);
N = N(1:150);
N = sort(N);

load('Training_test_data.mat')
y = [Data.y1(seed,N),Data.y2(seed,N)];
x = xreal(N);

yreal1 = Data.y1(seed,:);
yreal2 = Data.y2(seed,:);
y = y+randn(1,300)*0.1;

X = [x',ones(150,1); x',2*ones(150,1)];
t = X(:,2);
[Y,I] = sortrows(X,[2,1]);
t = t(I);
X = X(I,:);
Y = y';
Y = Y(I);

%Test locations
Xstar = [linspace(-8,8,3000)',ones(3000,1);linspace(-8,8,3000)',2*ones(3000,1)];

%Initial hyperparams
u = log(rand(6,1));
hyp.cov = [40000;0.5;0.5;-40000;0.5;0.5;u];
hyp.mean = mean(X(:,1));
hyp.lik = log(.2);

%Priors over hyperparams. Fix the branch/recomb. times. Everything else free.
p1 = {@priorClamped};
prior.mean = {[]};  % Gaussian prior for first, clamp second par
prior.lik = {[]};
prior.cov  = {p1;[];[];p1;[];[];[];[];[];[];[];[]};

im  = {@infPrior,@infExact,prior};                % inference method
par1 = {'meanConst','covBranchingRecombinationProcess_2A','likGauss',X,Y};
par2 = {'meanConst','covBranchingRecombinationProcess_2A','likGauss',X,Y,Xstar};
par3 = {'meanConst','covBranchingRecombinationProcess_2A','likGauss',X,Y};

hyp_pN      = feval(@minimize, hyp, @gp, -30000, im, par1{:});
[ymu ys2 fmu fs2   ] = gp(hyp_pN, im, par2{:});
[L dL   ] = gp(hyp_pN, im, par3{:});


if plt==1
 subplot(1,3,1);
 z = Xstar(1:3000,1);
 f = [fmu(1:3000)+2*sqrt(ys2(1:3000)); flipdim(fmu(1:3000)-2*sqrt(ys2(1:3000)),1)]; 
 fill([z; flipdim(z,1)], f, [7 7 7]/8);
 hold on
 z = Xstar(1:3000,1);
 f = [fmu(3000+1:end)+2*sqrt(ys2(3000+1:end)); flipdim(fmu(3000+1:end)-2*sqrt(ys2(3000+1:end)),1)]; 
 fill([z; flipdim(z,1)], f, [7 7 7]/8);
 plot(Xstar(1:3000,1),fmu(1:3000),'k-')
 plot(Xstar(3001:end,1),fmu(3001:end),'k-')
 plot(X(find(t==1),1),Y(find(t==1)),'ro'),
 plot(X(find(t==2),1),Y(find(t==2)),'bo')
end

Dist(1,:) = ((fmu(1:3000)-yreal1').^2);
Dist(2,:) = ((fmu(3001:end)-yreal2').^2);

%Individual fits ...
X1 = X(find(X(:,2)==1),1);
X2 = X(find(X(:,2)==2),1);
Y1 = Y(find(X(:,2)==1),1);
Y2 = Y(find(X(:,2)==2),1);
hyp.cov = u(1:2);%[log(rand(2,1))];%[3,0.5,0.5,1,3,0.3,2,2,2]';
hyp.mean = mean(X(find(X(:,2)==1),1));
hyp.lik = log(.2);


prior.mean = {[]};  % Gaussian prior for first, clamp second par
prior.lik = {[]};
prior.cov  = {[];[]};
im2  = {@infPrior,@infExact,prior}; 
par1 = {'meanConst','covSEiso','likGauss',X1,Y1};
par2 = {'meanConst','covSEiso','likGauss',X1,Y1,unique(Xstar(:,1))};
par3 = {'meanConst','covSEiso','likGauss',X1,Y1};

hyp_pN      = feval(@minimize, hyp, @gp, -30000, im2, par1{:});
[ymu1 ys21 fmu1 fs21   ] = gp(hyp_pN, im2, par2{:});
[L1 dL1   ] = gp(hyp_pN, im2, par3{:});

hyp.cov = u(3:4);%[3,0.5,0.5,1,3,0.3,2,2,2]';
hyp.mean = mean(X(find(X(:,2)==2),1));
hyp.lik = log(.2);

par1 = {'meanConst','covSEiso','likGauss',X2,Y2};
par2 = {'meanConst','covSEiso','likGauss',X2,Y2,unique(Xstar(:,1))};
par3 = {'meanConst','covSEiso','likGauss',X2,Y2};
hyp_pN      = feval(@minimize, hyp, @gp, -30000, im2, par1{:});
[ymu2 ys22 fmu2 fs22   ] = gp(hyp_pN, im2, par2{:});
[L2 dL2   ] = gp(hyp_pN, im2, par3{:});

Dist(3,:) = ((fmu1(1:3000)-yreal1').^2);
Dist(4,:) = ((fmu2(1:3000)-yreal2').^2);

%Joint GP
hyp.cov = u(1:2);%[3,0.5,0.5,1,3,0.3,2,2,2]';
hyp.mean = mean(X(find(X(:,2)==2),1));
hyp.lik = log(.2);
par1 = {'meanConst','covSEiso','likGauss',X(:,1),Y};
par2 = {'meanConst','covSEiso','likGauss',X(:,1),Y,unique(Xstar(:,1))};
par3 = {'meanConst','covSEiso','likGauss',X(:,1),Y};
hyp_pN      = feval(@minimize, hyp, @gp, -30000, im2, par1{:});
[ymu3 ys23 fmu3 fs23   ] = gp(hyp_pN, im2, par2{:});
[L3 dL3   ] = gp(hyp_pN, im2, par3{:});

Dist(5,:) = ((fmu1(1:3000)-yreal1').^2);
Dist(6,:) = ((fmu1(1:3000)-yreal2').^2);


Lout = [L,L1,L2,L3];


if plt==1
 subplot(1,3,2);
 z = Xstar(1:3000,1);
 f = [fmu1(1:3000)+2*sqrt(ys21(1:3000)); flipdim(fmu1(1:3000)-2*sqrt(ys21(1:3000)),1)]; 
 fill([z; flipdim(z,1)], f, [7 7 7]/8);
 hold on
 z = Xstar(1:3000,1);
 f = [fmu2(1:end)+2*sqrt(ys22(1:end)); flipdim(fmu2(1:end)-2*sqrt(ys22(1:end)),1)]; 
 fill([z; flipdim(z,1)], f, [7 7 7]/8);
 plot(Xstar(1:3000,1),fmu1(1:3000),'k-')
 plot(Xstar(1:3000,1),fmu2(1:end),'k-')
 plot(X(find(t==1),1),Y(find(t==1)),'ro'),
 plot(X(find(t==2),1),Y(find(t==2)),'bo')
 
  subplot(1,3,3);
 z = Xstar(1:3000,1);
 f = [fmu3(1:3000)+2*sqrt(ys23(1:3000)); flipdim(fmu3(1:3000)-2*sqrt(ys23(1:3000)),1)]; 
 fill([z; flipdim(z,1)], f, [7 7 7]/8);
 hold on
 plot(Xstar(1:3000,1),fmu1(1:3000),'k-')
 plot(Xstar(1:3000,1),fmu2(1:end),'k-')
 plot(X(find(t==1),1),Y(find(t==1)),'ro'),
 plot(X(find(t==2),1),Y(find(t==2)),'bo')
end


return