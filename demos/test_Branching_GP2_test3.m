function [Dist Lout] = test_Branching_GP2_test3(arc);

rng(arc)

addpath(genpath('/Users/christopher_penfold/Desktop/Code/gpss-research/source/gpml'))
addpath(genpath('/Users/christopher_penfold/Desktop/Code/deepGP/netlab3_3/'))

thetan = log(0.1);

xreal = linspace(-8,8,3000);
%xreal2 = [xreal',ones(3000,1);xreal',2*ones(3000,1)];
%K = feval('covHGP',log([0.5;2;0.1;2;1;2]'),xreal2);

%y = real(gsamp(zeros(1,6000),K,5));
%Data.y1 = y(:,1:3000);
%Data.y2 = y(:,3001:6000);
N = randperm(3000);
N = N(1:150);
N = sort(N);

load('HGPtraining.mat')
y = [Data.y1(arc,N),Data.y2(arc,N)];
x = xreal(N);

yreal1 = Data.y1(arc,:);
yreal2 = Data.y2(arc,:);

%yreal1 = zeros(1,3000);
%yreal1(find(xreal>=-pi/2 & xreal<=pi/2)) = cos(xreal(find(xreal>=-pi/2 & xreal<=pi/2)));
%yreal2 = zeros(1,3000);
%yreal2(find(xreal>=-pi/2 & xreal<=pi/2)) = -cos(xreal(find(xreal>=-pi/2 & xreal<=pi/2)));

% count = 0;
% x = randn(300000,1)*3;
% x = x(find(x<0));
% x = x(1:300);
% 
% for i = 1:length(x);
% if x(i)<-pi/2
%     y(i) = 0;
%     if rand(1,1)<0.5
%     t(i) = 1;
%     else
%     t(i) = 2;
%     end
% elseif x(i)>pi/2
% y(i) = 0;
% if rand(1,1)<0.5
% t(i) = 1;
% else
% t(i) = 2;
% end
% else
% 
% if rand(1)>0.5
% y(i) = cos(x(i));
% t(i) = 1;
% else
% y(i) = -cos(x(i));
% t(i) = 2;
% end
% end
% end
y = y+randn(1,300)*0.1;
%plot(x,y,'o')


%keyboard


%t = zeros(300,1);
%if  x(i)<-pi/2
%t(find(y>0))=2;
%t(find(y<=0))=1;

X = [x',ones(150,1); x',2*ones(150,1)];
t = X(:,2);
%[Y,I]=sort(X(:,2));
[Y,I] = sortrows(X,[2,1]);
t = t(I);
X = X(I,:);
%X = [x',ones(300,1);x',2*ones(300,1)];
%y = real(gsamp(zeros(1,300),Kall,10));
Y = y';
Y = Y(I);
%Y = y(1,1:300)' + rand(length(y(1,1:300)),1)*1;
Xstar = [linspace(-8,8,3000)',ones(3000,1);linspace(-8,8,3000)',2*ones(3000,1)];

u = log(rand(6,1));
hyp.cov = [4;0.5;0.5;-4;0.5;0.5;u];%[3,0.5,0.5,1,3,0.3,2,2,2]';
hyp.mean = mean(X(:,1));
hyp.lik = log(0.2);

hyp_pN      = feval(@minimize, hyp, @gp, -30000, @infExact, 'meanConst','covBranchingRecombinationProcess_v2','likGauss',X,Y);
[ymu ys2 fmu fs2   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingRecombinationProcess_v2', 'likGauss', X, Y, Xstar);
[L dL   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingRecombinationProcess_v2', 'likGauss', X, Y);
%[ymu ys2 fmu fs2   ] = gp(hyp, 'infExact', 'meanConst', 'covBranchingProcess', 'likGauss', X, Y, X);

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
hyp_pN      = feval(@minimize, hyp, @gp, -30000, @infExact, 'meanConst','covSEiso','likGauss',X1,Y1);
[ymu1 ys21 fmu1 fs21   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X1, Y1, unique(Xstar(:,1)));
[L1 dL1   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X1, Y1);

hyp.cov = u(1:2);%[3,0.5,0.5,1,3,0.3,2,2,2]';
hyp.mean = mean(X(find(X(:,2)==2),1));
hyp.lik = log(0.2);
hyp_pN      = feval(@minimize, hyp, @gp, -30000, @infExact, 'meanConst','covSEiso','likGauss',X2,Y2);
[ymu2 ys22 fmu2 fs22   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X2, Y2, unique(Xstar(:,1)));
[L2 dL2   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X2, Y2);

Dist(3,:) = ((fmu1(1:3000)-yreal1').^2);
Dist(4,:) = ((fmu2(1:3000)-yreal2').^2);



hyp.cov = u(1:2);%[3,0.5,0.5,1,3,0.3,2,2,2]';
hyp.mean = mean(X(:,1));
hyp.lik = log(0.2);
hyp_pN      = feval(@minimize, hyp, @gp, -30000, @infExact, 'meanConst','covSEiso','likGauss',X(:,1),Y);
[ymu3 ys23 fmu3 fs23   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X(:,1), Y(:,1), unique(Xstar(:,1)));
[L3 dL3   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covSEiso', 'likGauss', X(:,1), Y(:,1));

Dist(5,:) = ((fmu3-yreal1').^2);
Dist(6,:) = ((fmu3-yreal2').^2);


Lout = [L,L1,L2,L3];

%('HGP_accuracy_BRGP.mat')

% 
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
 f = [fmu3(1:3000)+2*sqrt(ys23(1:3000)); flipdim(fmu1(1:3000)-2*sqrt(ys23(1:3000)),1)]; 
 fill([z; flipdim(z,1)], f, [7 7 7]/8);
 hold on
 plot(Xstar(1:3000,1),fmu1(1:3000),'k-')
 plot(Xstar(1:3000,1),fmu2(1:end),'k-')
 plot(X(find(t==1),1),Y(find(t==1)),'ro'),
 plot(X(find(t==2),1),Y(find(t==2)),'bo')
%  


%Dist(3,:) = sqrt((fmu1(1:3000)-yreal1').^2);
%Dist(4,:) = sqrt((fmu2(1:3000)-yreal2').^2);

%keyboard
return
addpath('/Users/christopher_penfold/Downloads/ocsipackage/Functions/Other/')

%Now try to infer with wrong labels
Xtrue = X;
t = Xtrue(:,2);

X(find(X(:,1))>mean(X(:,1)),2) = 2;%randi([1 2],size(X,1),1);
X(find(X(:,1))<=mean(X(:,1)),2) = 1;

[Y,I]=sort(X(:,2));
X = X(I,:);
Y = Y(I);

hyp.cov = [1;0.5;0.5;-1;0.5;0.5;log(rand(6,1))];
hyp.mean = mean(X(:,1));
hyp.lik = 2;

hyp_pN      = feval(@minimize, hyp, @gp, -30000, @infExact, 'meanConst','covBranchingRecombinationProcess_v2','likGauss',X,Y);
%[ymu ys2 fmu fs2   ] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingRecombinationProcess_v2', 'likGauss', X, Y, Xstar);

[L dL] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingRecombinationProcess_v2', 'likGauss', X, Y);

Po(1,1) = L;
for i = 1:10000
    for j = 1:size(X,1)
        X1 = X;
        X2 = X;
        X1(j,2) = 1;
        X2(j,2) = 2;
    
        [Y11,I1]=sort(X1(:,2));
        [Y11,I2]=sort(X2(:,2));

        X1 = X1(I1,:);
        X2 = X2(I2,:);        
        Y1 = Y(I1);
        Y2 = Y(I2);
        
        [L1 dL] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingRecombinationProcess_v2', 'likGauss', X1, Y1);
        [L2 dL] = gp(hyp_pN, 'infExact', 'meanConst', 'covBranchingRecombinationProcess_v2', 'likGauss', X2, Y2);

        po = [L1,L2];
        po = exp(-po)./sum(exp(-po));
        s = discreternd(po,1);

        if s==1
            X = X1;
            Y = Y1;
            Po(i,1) = L1;
        else
            X = X2;
            Y = Y2;
            Po(i,1) = L2;            
        end
        
    end    
end