%Demo plot showing different examples of branching/recombination GPs possible 
%via composition of different base kernels with changepoint kernels.

addpath(genpath('../'))
%addpath(genpath('/Users/christopher_penfold/Desktop/Code/gpss-research/source/gpml'))
%addpath(genpath('/Users/christopher_penfold/Desktop/Code/deepGP/netlab3_3/'))
x = linspace(0,10,100);

covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
K1   = covSEiso([1,3],x');
K2   = feval(covN1{:}, [3,0.5,0.2,3], x');
Kall = [K1,K1;K1,K1+K2];
y=real(gsamp(zeros(1,200),Kall,10));

hFig1 = figure(1)
subplot(4,3,3);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,real(y(1,201:300))','g')
subplot(4,3,2);plot(x,real(y(2,1:100))','b'),hold on, plot(x,real(y(2,101:200))','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,real(y(2,201:300))','g')
%title('A smooth function branches from another smooth function')
subplot(4,3,1); imagesc(Kall),,set(gca,'XTick',[],'YTick',[])


covN1 = {@covChangePointMultiD, {1, @covZero, @covPeriodic}};
K1   = covSEiso([1,3],x');
K2   = feval(covN1{:}, [3,0.5,0.5,0.2,3], x');
%K3   = feval(covN1{:}, [3,0.5,0.5,0.6,3], x');
Kall = [K1,K1;K1,K1+K2];

y=real(gsamp(zeros(1,200),Kall,10));
subplot(4,3,6);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,real(y(1,201:300))','g')
subplot(4,3,5);plot(x,real(y(2,1:100))','b'),hold on, plot(x,real(y(2,101:200))','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,real(y(2,201:300))','g')
%title('Oscillating function branches from smooth function')
subplot(4,3,4); imagesc(Kall),set(gca,'XTick',[],'YTick',[])

covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
%covN1 = {@covChangePointMultiD, {1, @covZero, @covNNoone}};
K1   = covSEiso([1,3],x');
K2   = feval(covN1{:}, [3,0.5,1,2], x');
K3   = feval(covN1{:}, [7,0.5,3,3], x');
Kall = [K1,K1,K1;
        K1,K1+K2,K1+K2;
        K1,K1+K2,K1+K2+K3];
y=real(gsamp(zeros(1,300),Kall,10));


subplot(4,3,12);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),,plot(x,real(y(1,201:300))','k-'),set(gca,'XTick',[],'YTick',[])
subplot(4,3,11);plot(x,real(y(2,1:100))','b'),hold on, plot(x,real(y(2,101:200))','r'),plot(x,real(y(2,201:300))','k-'),set(gca,'XTick',[],'YTick',[])
%title('Two smooth functions branched from a latent smooth function')
subplot(4,3,10); imagesc(Kall),set(gca,'XTick',[],'YTick',[])



covN1 = {@covChangePointMultiD, {1, @covZero, @covPeriodic}};
covN2 = {@covChangePointMultiD, {1, @covZero, @covPeriodic}};

K1   = covSEiso([1,3],x');
K2   = feval(covN1{:}, [3,0.5,0.5,0.2,3], x');
K3   = feval(covN1{:}, [7,0.5,0.5,0.6,3], x');
Kall = [K1,K1,K1;K1,K1+K2,K1;K1,K1,K1+K3];

y=real(gsamp(zeros(1,300),Kall,10));
subplot(4,3,9);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),plot(x,real(y(1,201:300))','k-'),set(gca,'XTick',[],'YTick',[])
subplot(4,3,8);plot(x,real(y(2,1:100))','b'),hold on, plot(x,real(y(2,101:200))','r'),plot(x,real(y(2,201:300))','k-'),set(gca,'XTick',[],'YTick',[])
%title('Two oscillating functions branched from latent smooth function')
subplot(4,3,7); imagesc(Kall),set(gca,'XTick',[],'YTick',[])
hFig1.PaperUnits = 'centimeters';
hFig1.PaperPosition = [0 0 8.6 15.0];
print('./results/plots/Plot1','-dpng','-r0')


%covN2 = {@covChangePointMultiD, {1, {@covChangePointMultiD, {1, @covZero, @covNoise}}, @covZero}};

% %covN1 = {@covChangePointMultiD,{1,@covZero,{@covChangePointMultiD, {1, @covNoise, @covZero}}};
% %covN2 = {@covChangePointMultiD,{1,@covZero,{@covChangePointMultiD, {1, @covNoise, @covZero}}};
% 
% K1   = covSEiso([1,3],x');
% K2   = feval(covN2{:}, [6,0.5,4,0.2,2.5], x');
% K3   = feval(covN2{:}, [6,0.5,4,0.2,1], x');
% %K3   = feval(covN1{:}, [3,0.5,0.5,0.6,3], x');
% Kall = [K2+K1,K1,K1;K1,K1+K3,K1;K1,K1,K1];
% 
% y=real(gsamp(zeros(1,300),Kall,10));
% 
% subplot(5,3,15);plot(x,y(1,1:100)','b'),hold on, plot(x,y(1,101:200)','r'),set(gca,'XTick',[],'YTick',[]),set(gca,'XTick',[],'YTick',[])%,plot(x,y(1,201:300)','k'),plot(x,y(1,301:400)','g')
% subplot(5,3,14);plot(x,y(3,1:100)','b'),hold on, plot(x,y(3,101:200)','r'),set(gca,'XTick',[],'YTick',[]),set(gca,'XTick',[],'YTick',[])%,plot(x,y(3,201:300)','k'),plot(x,y(3,301:400)','g')
% title('Identical functions with varied locus-specific white noise')
% subplot(5,3,13); imagesc(Kall),set(gca,'XTick',[],'YTick',[])
% 



hFig1 =figure(2)
covN1 = {@covChangePointMultiD, {1, @covSEiso,@covZero}};
K1   = covSEiso([1,3],x');
K2   = feval(covN1{:}, [3,0.5,0.2,3], x');
Kall = [K1,K1;K1,K1+K2];
y=real(gsamp(zeros(1,200),Kall,10));


subplot(4,3,3);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,real(y(1,201:300))','g')
subplot(4,3,2);plot(x,real(y(2,1:100))','b'),hold on, plot(x,real(y(2,101:200))','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,real(y(2,201:300))','g')
%title('A smooth function branches from another smooth function')
subplot(4,3,1); imagesc(Kall),,set(gca,'XTick',[],'YTick',[])


covN1 = {@covChangePointMultiD, {1, @covPeriodic, @covZero}};
K1   = covSEiso([1,3],x');
K2   = feval(covN1{:}, [3,0.5,0.5,0.2,3], x');
%K3   = feval(covN1{:}, [3,0.5,0.5,0.6,3], x');
Kall = [K1,K1;K1,K1+K2];

y=real(gsamp(zeros(1,200),Kall,10));
subplot(4,3,6);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,real(y(1,201:300))','g')
subplot(4,3,5);plot(x,real(y(2,1:100))','b'),hold on, plot(x,real(y(2,101:200))','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,real(y(2,201:300))','g')
%title('Oscillating function branches from smooth function')
subplot(4,3,4); imagesc(Kall),set(gca,'XTick',[],'YTick',[])


covN1 = {@covChangePointMultiD, {1, @covSEiso, @covZero}};
K1   = covSEiso([1,3],x');
K2   = feval(covN1{:}, [7,0.5,1,2], x');
K3   = feval(covN1{:}, [3,0.5,3,3], x');
Kall = [K1,K1,K1;K1,K1+K2,K1;K1,K1,K1+K3];
y=real(gsamp(zeros(1,300),Kall,10));


subplot(4,3,9);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),,plot(x,real(y(1,201:300))','k-'),set(gca,'XTick',[],'YTick',[])
subplot(4,3,8);plot(x,real(y(2,1:100))','b'),hold on, plot(x,real(y(2,101:200))','r'),plot(x,real(y(2,201:300))','k-'),set(gca,'XTick',[],'YTick',[])
%title('Two smooth functions branched from a latent smooth function')
subplot(4,3,7); imagesc(Kall),,set(gca,'XTick',[],'YTick',[])



covN1 = {@covChangePointMultiD, {1, @covPeriodic,@covZero}};
covN2 = {@covChangePointMultiD, {1, @covPeriodic,@covZero}};

K1   = covSEiso([1,3],x');
K2   = feval(covN1{:}, [7,0.5,0.5,0.2,3], x');
K3   = feval(covN1{:}, [3,0.5,0.5,0.6,3], x');
Kall = [K1,K1,   K1;
        K1,K1+K2,K1+K2;
        K1,K1+K2,K1+K2+K3];

y=real(gsamp(zeros(1,300),Kall,10));
subplot(4,3,12);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),plot(x,real(y(1,201:300))','k-'),set(gca,'XTick',[],'YTick',[])
subplot(4,3,11);plot(x,real(y(2,1:100))','b'),hold on, plot(x,real(y(2,101:200))','r'),plot(x,real(y(2,201:300))','k-'),set(gca,'XTick',[],'YTick',[])
%title('Two oscillating functions branched from latent smooth function')
subplot(4,3,10); imagesc(Kall),set(gca,'XTick',[],'YTick',[])


hFig1.PaperUnits = 'centimeters';
hFig1.PaperPosition = [0 0 8.6 15.0];
print('./results/plots/Plot2','-dpng','-r0')


hFig1 = figure(3)

covN2 = {@covChangePointMultiD, {1, {@covChangePointMultiD, {1, @covZero, @covPeriodic}}, @covZero}};

K2   = feval(covN2{:}, [8,0.5,3,0.5,0.5,0.2,3], x');
Kall = [K2+K1,K1;K1,K1];

y=real(gsamp(zeros(1,200),Kall,10));
subplot(4,3,3);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),set(gca,'XTick',[],'YTick',[])
subplot(4,3,2);plot(x,real(y(2,1:100))','b'),hold on, plot(x,real(y(2,101:200))','r'),set(gca,'XTick',[],'YTick',[])
%title('Oscillating function branches and recombines with smooth function')
subplot(4,3,1); imagesc(Kall),set(gca,'XTick',[],'YTick',[])


K1a   = covPeriodic([0.5,0.3,3],x');
K2   = feval(covN2{:}, [8,0.5,3,0.5,0.5,0.2,3], x');
Kall = [K2+K1a,K1a;K1a,K1a];

y=real(gsamp(zeros(1,200),Kall,10));
subplot(4,3,6);plot(x,real(y(1,1:100))','b'),hold on, plot(x,real(y(1,101:200))','r'),set(gca,'XTick',[],'YTick',[])
subplot(4,3,5);plot(x,real(y(2,1:100))','b'),hold on, plot(x,real(y(2,101:200))','r'),set(gca,'XTick',[],'YTick',[])
%title('Two oscillating function branch and recombine')
subplot(4,3,4); imagesc(Kall),set(gca,'XTick',[],'YTick',[])

% 
 covN2 = {@covChangePointMultiD, {1, {@covChangePointMultiD, {1, @covZero, @covSEiso}}, @covZero}};
% K2   = feval(covN2{:}, [8,0.5,3,0.5,1,3], x');
% K3   = feval(covN2{:}, [6,0.5,5.5,0.5,1,3], x');
% 
% Kall = [K1,K1,K1;K1,K1+K2,K1+K2;K1,K1+K2,K1+K2+K3];
% 
% y=real(gsamp(zeros(1,300),Kall,10));
% subplot(5,3,9);plot(x,y(1,1:100)','b'),hold on, plot(x,y(1,101:200)','r'),plot(x,y(1,201:300)','k')
% subplot(5,3,8);plot(x,y(3,1:100)','b'),hold on, plot(x,y(3,101:200)','r'),plot(x,y(3,201:300)','k')
% %title('Set of smooth functions branch and recombines with smooth function')
% subplot(5,3,7); imagesc(Kall)

K2   = feval(covN2{:}, [8,0.5,3,0.5,1,3], x');
K3   = feval(covN2{:}, [6,0.5,5.5,0.5,1,3], x');
K4   = feval(covN2{:}, [4,0.5,1,0.5,1,3], x');
Kall = [K1,K1,   K1,      K1;
        K1,K1+K2,K1+K2,   K1;
        K1,K1+K2,K1+K2+K3,K1;
        K1,K1   ,K1,      K1+K4];
        
y=real(gsamp(zeros(1,400),Kall,10));
subplot(4,3,9);plot(x,y(1,1:100)','b'),hold on, plot(x,y(1,101:200)','r'),plot(x,y(1,201:300)','k'),plot(x,y(1,301:400)','g'),set(gca,'XTick',[],'YTick',[])
subplot(4,3,8);plot(x,y(3,1:100)','b'),hold on, plot(x,y(3,101:200)','r'),plot(x,y(3,201:300)','k'),plot(x,y(3,301:400)','g'),set(gca,'XTick',[],'YTick',[])
%title('Set of smooth functions branch and recombines with smooth function')
subplot(4,3,7); imagesc(Kall),set(gca,'XTick',[],'YTick',[])




covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
covN2 = {@covChangePointMultiD, {1, @covSEiso, @covZero}};
covN3 = {@covChangePointMultiD, {1, {@covChangePointMultiD, {1, @covZero, @covSEiso}}, @covZero}};

K2   = feval(covN1{:}, [4,0.5,1,3], x');
K3   = feval(covN2{:}, [5,0.5,1,3], x');
K4   = feval(covN3{:}, [6,0.5,2,0.5,1,3], x');

Kall = [K1,K1,   K1,      K1;
        K1,K1+K2,K1+K2,   K1+K2;
        K1,K1+K2,K1+K2+K3,K1+K2+K3;
        K1,K1+K2   ,K1+K2+K3,      K1+K2+K3+K4];
y=real(gsamp(zeros(1,400),Kall,10));
subplot(4,3,12);plot(x,y(1,1:100)','b'),hold on, plot(x,y(1,101:200)','r'),plot(x,y(1,201:300)','k'),plot(x,y(1,301:400)','g'),set(gca,'XTick',[],'YTick',[])
subplot(4,3,11);plot(x,y(3,1:100)','b'),hold on, plot(x,y(3,101:200)','r'),plot(x,y(3,201:300)','k'),plot(x,y(3,301:400)','g'),set(gca,'XTick',[],'YTick',[])
%title('Set of smooth functions branch and recombines with smooth function')
subplot(4,3,10); imagesc(Kall),set(gca,'XTick',[],'YTick',[])


hFig1.PaperUnits = 'centimeters';
hFig1.PaperPosition = [0 0 8.6 15.0];
print('./results/plots/Plot3','-dpng','-r0')

hFig1 = figure(4)
covN2 = {@covChangePointMultiD, {1, {@covChangePointMultiD, {1, @covZero, @covNoise}}, @covZero}};
K1   = covSEiso([1,3],x');
K2   = feval(covN2{:}, [6,0.5,4,0.2,2.5], x');
K3   = feval(covN2{:}, [6,0.5,4,0.2,1], x');
Kall = [K2+K1,K1,K1;K1,K1+K3,K1;K1,K1,K1];

y=real(gsamp(zeros(1,300),Kall,10));
subplot(1,3,3);plot(x,y(1,1:100)','b'),hold on, plot(x,y(1,101:200)','r'),set(gca,'XTick',[],'YTick',[]),set(gca,'XTick',[],'YTick',[])%,plot(x,y(1,201:300)','k'),plot(x,y(1,301:400)','g')
subplot(1,3,2);plot(x,y(3,1:100)','b'),hold on, plot(x,y(3,101:200)','r'),set(gca,'XTick',[],'YTick',[]),set(gca,'XTick',[],'YTick',[])%,plot(x,y(3,201:300)','k'),plot(x,y(3,301:400)','g')
%title('Identical functions with varied locus-specific white noise')
subplot(1,3,1); imagesc(Kall),set(gca,'XTick',[],'YTick',[])


% covN1 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
% covN2 = {@covChangePointMultiD, {1, @covSEiso, @covZero}};
% covN3 = {@covChangePointMultiD, {1, {@covChangePointMultiD, {1, @covZero, @covSEiso}}, @covZero}};
% 
% K2   = feval(covN1{:}, [4,0.5,1,3], x');
% K3   = feval(covN2{:}, [5,0.5,1,3], x');
% K4   = feval(covN3{:}, [6,0.5,2,0.5,1,3], x');
% 
% Kall = [K1,K1,   K1,      K1;
%         K1,K1+K2,K1+K2,   K1+K2;
%         K1,K1+K2,K1+K2+K3,K1+K2+K3;
%         K1,K1+K2   ,K1+K2+K3,      K1+K2+K3+K4];
% y=real(gsamp(zeros(1,400),Kall,10));
% subplot(2,3,6);plot(x,y(1,1:100)','b'),hold on, plot(x,y(1,101:200)','r'),plot(x,y(1,201:300)','k'),plot(x,y(1,301:400)','g')
% subplot(2,3,5);plot(x,y(3,1:100)','b'),hold on, plot(x,y(3,101:200)','r'),plot(x,y(3,201:300)','k'),plot(x,y(3,301:400)','g')
% %title('Set of smooth functions branch and recombines with smooth function')
% subplot(2,3,4); imagesc(Kall)

hFig1.PaperUnits = 'centimeters';
hFig1.PaperPosition = [0 0 8.6 3.0];
print('./results/plots/Plot4','-dpng','-r0')

return

covN2 = {@covChangePointMultiD, {1, {@covChangePointMultiD, {1, @covZero, @covNoise}}, @covZero}};

%covN1 = {@covChangePointMultiD,{1,@covZero,{@covChangePointMultiD, {1, @covNoise, @covZero}}};
%covN2 = {@covChangePointMultiD,{1,@covZero,{@covChangePointMultiD, {1, @covNoise, @covZero}}};

K1   = covSEiso([1,3],x');
K2   = feval(covN2{:}, [6,0.5,4,0.2,2.5], x');
K3   = feval(covN2{:}, [6,0.5,4,0.2,1], x');
%K3   = feval(covN1{:}, [3,0.5,0.5,0.6,3], x');
Kall = [K2+K1,K1,K1;K1,K1+K3,K1;K1,K1,K1];

y=real(gsamp(zeros(1,300),Kall,10));

figure(6)
subplot(1,3,3);plot(x,y(1,1:100)','b'),hold on, plot(x,y(1,101:200)','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,y(1,201:300)','k'),plot(x,y(1,301:400)','g')
subplot(1,3,2);plot(x,y(3,1:100)','b'),hold on, plot(x,y(3,101:200)','r'),set(gca,'XTick',[],'YTick',[])%,plot(x,y(3,201:300)','k'),plot(x,y(3,301:400)','g')
title('Set of smooth functions branch and recombines with smooth function')
subplot(1,3,1); imagesc(Kall)
