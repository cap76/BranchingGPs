%Demo plots using Branching GPs for transitions models using branching of constant kernels. 
%Potential uses include change-point/switch models of gene expression in which gene expression
%transitions between steady states

addpath(genpath('../'))
%addpath(genpath('/Users/christopher_penfold/Desktop/Code/gpss-research/source/gpml'))
%addpath(genpath('/Users/christopher_penfold/Desktop/Code/deepGP/netlab3_3/'))

x = linspace(0,10,100);

%Define a branching process between constant covariances
covN1 = {@covChangePointMultiD, {1, @covZero, @covConst}};
K1   = covConst(0.1,x');
K2   = feval(covN1{:}, [3,2,2], x');
Kall = [K2+K1,K1;K1,K1];
y=real(gsamp(zeros(1,200),Kall,10));

hFig1 = figure(1)
subplot(2,2,2);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,real(y(1,201:300))','g')
subplot(2,2,1); imagesc(Kall),set(gca,'XTick',[],'YTick',[])

%Define transitions in 2D space
x = linspace(0,10,50);
y = linspace(0,10,50);
X = repmat(x,50,1);
Y = repmat(y',1,50);
Xin = [reshape(X,50^2,1),reshape(Y,50^2,1)];

K1 = covConst([1],Xin);
covN1 = {@covProd,{{@covChangePointMultiD, {1, @covZero, @covConst}},{@covChangePointMultiD, {2, @covZero, @covConst}}}};

K2   = feval(covN1{:}, [3,0.5,2,3,0.5,2], Xin);
Kall = [K2+K1,K1;K1,K1];

y=real(gsamp(zeros(1,2* 2500),Kall,10));
subplot(2,2,4);plot3(Xin(:,1),Xin(:,2),real(y(2,1:2500))','r.')
hold on
plot3(Xin(:,1),Xin(:,2),real(y(2,2501:end))','b.')

hFig1.PaperUnits = 'centimeters';
hFig1.PaperPosition = [0 0 8.6 3.0];
print('./results/plots/Plot5','-dpng','-r0')

