function [Dist Lout] = test_Branching_GP2_test1(seed,plt);
%Benchmark a branching/recombination GP and compare ML/SSE to independent GPs and joint GPs.
%Input: seed number (seed) and plot (plt). To benchmark loop seeds 1-100. 

rng(seed)

addpath(genpath('../'))
%addpath(genpath('/Users/christopher_penfold/Desktop/Code/gpss-research/source/gpml'))
%addpath(genpath('/Users/christopher_penfold/Desktop/Code/deepGP/netlab3_3/'))

%Generate data from a simple branching/recombination process on the fly.
xreal = linspace(-8,8,3000);
yreal1 = zeros(1,3000);
yreal1(find(xreal>=-pi/2 & xreal<=pi/2)) = cos(xreal(find(xreal>=-pi/2 & xreal<=pi/2)));
yreal2 = zeros(1,3000);
yreal2(find(xreal>=-pi/2 & xreal<=pi/2)) = -cos(xreal(find(xreal>=-pi/2 & xreal<=pi/2)));
count = 0;
x = randn(300,1)*3;
for i = 1:length(x);
if x(i)<-pi/2
    y(i) = 0;
    if rand(1,1)<0.5
    t(i) = 1;
    else
    t(i) = 2;
    end
elseif x(i)>pi/2
y(i) = 0;
if rand(1,1)<0.5
t(i) = 1;
else
t(i) = 2;
end
else

if rand(1)>0.5
y(i) = cos(x(i));
t(i) = 1;
else
y(i) = -cos(x(i));
t(i) = 2;
end
end
end
y = y+randn(1,300)*0.1;

X = [x,t'];

%Order data
[Y,I] = sortrows(X,[2,1]);
t = t(I);
X = X(I,:);
Y = y';
Y = Y(I);

%Prediciton points
Xstar = [linspace(-8,8,3000)',ones(3000,1);linspace(-8,8,3000)',2*ones(3000,1)];

%Initial hyperparameters (these will be shared across models to make test fairer)
u = log(rand(6,1));
hyp.cov = [4;0.5;0.5;-4;0.5;0.5;u];
hyp.mean = mean(X(:,1));
hyp.lik = log(.2);

%Branching/recombination process opt.
hyp_pN      = feval(@minimize, hyp, @gp, -30000, @infExact, 'meanConst','covBranchingRecombinationProcess_2A','likGauss',X,Y);
[ymu ys2 fmu fs2   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingRecombinationProcess_2A', 'likGauss', X, Y, Xstar);
[L dL   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingRecombinationProcess_2A', 'likGauss', X, Y);

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

%Individual GP fits ...
X1 = X(find(X(:,2)==1),1);
X2 = X(find(X(:,2)==2),1);
Y1 = Y(find(X(:,2)==1),1);
Y2 = Y(find(X(:,2)==2),1);
hyp.cov = u(1:2);%[log(rand(2,1))];%[3,0.5,0.5,1,3,0.3,2,2,2]';
hyp.mean = mean(X(find(X(:,2)==1),1));
hyp.lik = log(.2);
hyp_pN      = feval(@minimize, hyp, @gp, -30000, @infExact, 'meanConst','covSEiso','likGauss',X1,Y1);
[ymu1 ys21 fmu1 fs21   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X1, Y1, unique(Xstar(:,1)));
[L1 dL1   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X1, Y1);

hyp.cov = u(3:4);
hyp.mean = mean(X(find(X(:,2)==2),1));
hyp.lik = log(.2);
hyp_pN      = feval(@minimize, hyp, @gp, -30000, @infExact, 'meanConst','covSEiso','likGauss',X2,Y2);
[ymu2 ys22 fmu2 fs22   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X2, Y2, unique(Xstar(:,1)));
[L2 dL2   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X2, Y2);

Dist(3,:) = ((fmu1(1:3000)-yreal1').^2);
Dist(4,:) = ((fmu2(1:3000)-yreal2').^2);

%Joint GP fit
hyp.cov = u(1:2);
hyp.mean = mean(X(:,1));
hyp.lik = log(.2);
hyp_pN      = feval(@minimize, hyp, @gp, -30000, @infExact, 'meanConst','covSEiso','likGauss',X(:,1),Y);
[ymu3 ys23 fmu3 fs23   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X(:,1), Y, unique(Xstar(:,1)));
[L3 dL3   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X(:,1), Y);

Dist(5,:) = ((fmu3(1:3000)-yreal1').^2);
Dist(6,:) = ((fmu3(1:3000)-yreal2').^2);

%Return marginal likeihoods
Lout = [L,L1,L2,L3];

if plt==1
 %Plot fits
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