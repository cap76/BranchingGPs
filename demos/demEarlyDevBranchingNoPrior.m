function [Fin] = demEarlyDevBranchingNoPrior(seed,batchi);


warning off all
%addpath(genpath('../Pseudotime_Algorithms/Functions/gpml-matlab-v3.6-2015-07-07/'))
%addpath(('./gpml/cov'))
%addpath(genpath('./netlab3_3/'))
addpath(genpath('../'))

rng(seed) %Ensure the 

%load('PGClabelleddata.mat');
D1 = importdata('./AllProcesses_earlydev.csv');
%D1 = importdata('/Users/christopher_penfold/Desktop/Pathway Modelling/Datasets/AllDE_P_gtp1_FC_gt10.csv')

Type = D1.data(1,:);%D1.data(1,:);
CapT = D1.data(2,:);
Time = D1.data(3,:);%D1.data(3,:);
%ind = find(Time>-1);
Data.Y = log2(D1.data(4:end,:)+1)';
%Time = Time(ind);
%Sex = Sex(ind);
%Type = Type(ind);

z = zeros(1,length(CapT));
z(find(Type==2))=ones;
z(find(Type==3))=2*ones;
z(find(Type==1))=3*ones;
%z(find(Type==1))=4*ones;
z(find(z==0))=randi([1 3],1,length(find(z==0)));

remove_inds = [];
uz = unique(Type);
ut = unique(CapT);
for i = 1:length(ut)
    for j = 1:length(uz)
        ind = find(CapT==ut(i) & Type==uz(j));
        I(i,j) = length(ind);

        if isempty(ind)==0                        
                    MT(i,j) = mean(Time(ind));
        remove_ind = ind(randperm(length(ind)));
        remove_inds = [remove_inds,remove_ind(1:length(ind)-20)];
        end
        
    end
end

ind_train = setdiff(1:1:length(Time),[remove_inds]);

Type_pr = Type(remove_inds);
CapT_pr = CapT(remove_inds);
Time_pr = Time(remove_inds);
Y_pr    = Data.Y(remove_inds,:);
z_pr    = z(remove_inds);

Type    = Type(ind_train);
CapT    = CapT(ind_train);
Time    = Time(ind_train);
Data.Y  = Data.Y(ind_train,:);
z       = z(ind_train);

Data.X = [Time',z'];

%1,2,3 -> 1,2,2
Data.X1 = Data.X;
Data.X1(find(Data.X1(:,2)==3),2)=2;

%1,2,3 -> 1,1,2
Data.X2 = Data.X;
Data.X2(find(Data.X2(:,2)==2),2)=1;
Data.X2(find(Data.X2(:,2)==3),2)=2;

%1,2,3 -> 1,2,1
Data.X3 = Data.X;
Data.X3(find(Data.X3(:,2)==3),2)=1;

%1,2,3 -> 1,1,1
Data.X4 = Data.X;
Data.X4(find(Data.X4(:,2)==2),2)=1;
Data.X4(find(Data.X4(:,2)==3),2)=1;
%for i = 1:size(Data.Y,2)
%    Data.Y(:,i) = (Data.Y(:,i)-mean(Data.Y(:,i)))
%end
%keyboard

%D1 = importdata('./Pa13-Combo2.txt');

try %Try to resume analysis 
    load(['../results/earlydevelopment/EarlDevResults_' num2str(batchi) '_' num2str(seed) '_noprior_r.mat'])
    startind = (batchi-1)*100 + 1;%startind = length(Ps)+1;
    endind   = (batchi)*100;
catch
    startind = (batchi-1)*100 + 1;
    endind   = (batchi)*100;
end

endind = min(endind,size(Data.Y,2));

%X1 = [repmat(1:1:13,1,12); ones(1,52),2*ones(1,52),3*ones(1,52)]';
%X2 = [repmat(1:1:13,1,12); ones(1,52),ones(1,52),2*ones(1,52)]';
%X3 = [repmat(1:1:13,1,12); ones(1,52),2*ones(1,52),2*ones(1,52)]';

Xstar1 = [repmat(linspace(0,30,1000),1,4);ones(1,1000),2*ones(1,1000),3*ones(1,1000),4*ones(1,1000)]';

%Xstar2 = [repmat(linspace(0,30,1000),1,4);ones(1,1000),ones(1,1000),2*ones(1,1000),4*ones(1,1000)]';
%Xstar3 = [repmat(linspace(0,30,1000),1,4);ones(1,1000),2*ones(1,1000),2*ones(1,1000),4*ones(1,1000)]';
%Xstar4 = [repmat(linspace(0,30,1000),1,4);2*ones(1,1000),ones(1,1000),2*ones(1,1000),4*ones(1,1000)]';
%Xstar5 = [repmat(linspace(0,30,1000),1,4);ones(1,1000),ones(1,1000),ones(1,1000),4*ones(1,1000)]';

%Xstar2 = [repmat(linspace(0,13,50),1,3);ones(1,50),1*ones(1,50),2*ones(1,50)]';
%Xstar3 = [repmat(linspace(0,13,50),1,3);ones(1,50),2*ones(1,50),2*ones(1,50)]';

%keyboard

for i = startind:endind

    %try
    %[length(Ps),startind]
    %catch
    %end
    
%Y1 = D1.data(i,:)'; %Mock hrp DC
%Y2 = D1.data(i,[1:52,2*52+1:3*52,52+1:2*52])';
%l1 = log(3); l2 = log(3); lg = log(3); v1 = log(3); v2 = log(3); vg = log(std(Data.Y(:,i)));

%1->2->3. Mock->hrp->DC
%keyboard
%
%u = linspace(min(Data.X(:,1)),max(Data.X(:,1)),30);
%uhat = [u',ones(30,1);u',2*ones(30,1);u',3*ones(30,1)];
%covfuncF = {@covFITC, {@covBranchingProcess_4B}, uhat};




%hyp.xu = uhat;
hyp.cov  = [11;1.5;11;1.5;11;2;2;2;2;2;2;2;2;2]; hyp.mean = mean(Data.Y(:,i)); hyp.lik  = 2;

%Use : covBranchingProcess_3C

 hyp_pN1 = feval(@minimize, hyp, @gp, -20000, 'infExact', 'meanConst', @covBranchingProcess_3C,'likGauss',Data.X,Data.Y(:,i));         % optimise
 [L1 dL1] = feval(@gp,hyp_pN1, 'infExact', 'meanConst',@covBranchingProcess_3C,'likGauss',Data.X,Data.Y(:,i));         % optimise
 [ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, 'infExact', 'meanConst', @covBranchingProcess_3C,'likGauss',Data.X,Data.Y(:,i),Xstar1);
 
 
 hyp_pN2 = feval(@minimize, hyp, @gp, -20000, 'infExact', 'meanConst', @covBranchingProcess_3C,'likGauss',Data.X1,Data.Y(:,i));         % optimise
 [L2 dL2] = feval(@gp,hyp_pN2, 'infExact', 'meanConst',@covBranchingProcess_3C,'likGauss',Data.X1,Data.Y(:,i));         % optimise
 [ymu2 ys22 fmu2 fs22   ]= feval(@gp,hyp_pN2, 'infExact', 'meanConst', @covBranchingProcess_3C,'likGauss',Data.X1,Data.Y(:,i),Xstar1);
 
 hyp_pN3 = feval(@minimize, hyp, @gp, -20000, 'infExact', 'meanConst', @covBranchingProcess_3C,'likGauss',Data.X2,Data.Y(:,i));         % optimise
 [L3 dL3] = feval(@gp,hyp_pN3, 'infExact', 'meanConst',@covBranchingProcess_3C,'likGauss',Data.X2,Data.Y(:,i));         % optimise
 [ymu3 ys23 fmu3 fs23   ]= feval(@gp,hyp_pN3, 'infExact', 'meanConst', @covBranchingProcess_3C,'likGauss',Data.X2,Data.Y(:,i),Xstar1);
 
 
 hyp_pN4 = feval(@minimize, hyp, @gp, -20000, 'infExact', 'meanConst', @covBranchingProcess_3C,'likGauss',Data.X3,Data.Y(:,i));         % optimise
 [L4 dL4] = feval(@gp,hyp_pN4, 'infExact', 'meanConst',@covBranchingProcess_3C,'likGauss',Data.X3,Data.Y(:,i));         % optimise
 [ymu4 ys24 fmu4 fs24   ]= feval(@gp,hyp_pN4, 'infExact', 'meanConst', @covBranchingProcess_3C,'likGauss',Data.X3,Data.Y(:,i),Xstar1);
 


 %hyp_pN5 = feval(@minimize, hyp, @gp, -20000, 'infExact', 'meanConst', @covBranchingProcess_3C,'likGauss',Data.X4,Data.Y(:,i));         % optimise
 %[L5 dL5] = feval(@gp,hyp_pN5, 'infExact', 'meanConst',@covBranchingProcess_3C,'likGauss',Data.X4,Data.Y(:,i));         % optimise
 %[ymu5 ys25 fmu5 fs25   ]= feval(@gp,hyp_pN5, 'infExact', 'meanConst', @covBranchingProcess_3C,'likGauss',Data.X4,Data.Y(:,i),Xstar1);


%hyp_pN1 = feval(@minimize, hyp, @gp, -10000, 'infFITC', 'meanConst', covfuncF,'likGauss',Data.X,Data.Y(:,i));         % optimise
%[L1 dL1] = feval(@gp,hyp_pN1, 'infFITC', 'meanConst',covfuncF,'likGauss',Data.X,Data.Y(:,i));         % optimise
%[ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, 'infExact', 'meanConst','covBranchingProcess_v4c','likGauss',Data.X,Data.Y(:,i),Xstar1);

%1->2,1->3, Mock->hrp, Mock->DC
%hyp_pN2 = feval(@minimize, hyp, @gp, -10000, 'infExact','meanConst','covBranchingProcess_v4','likGauss',Data.X1,Data.Y1(:,i));         % optimise
%[L2 dL2] = feval(@gp,hyp_pN2,'infExact','meanConst','covBranchingProcess_v4','likGauss',Data.X1,Data.Y1(:,i));         % optimise
%[ymu2 ys22 fmu2 fs22   ]= feval(@gp,hyp_pN2, 'infExact','meanConst','covBranchingProcess_v4','likGauss',Data.X1,Data.Y1(:,i),Xstar1);
%keyboard
%Now do joint modelling of subsets
%1,2->3, Mock,Hrp -> DC
%hyp.cov = [2;0.5;l1;v1;lg;vg];hyp.mean = mean(Data.Y1(:,i));hyp.lik = 2;
%hyp_pN3 = feval(@minimize, hyp, @gp, -10000, 'infExact','meanConst','covBranchingProcess_v5','likGauss',Data.X3,Data.Y3(:,i));         % optimise
%[L3 dL3] = feval(@gp,hyp_pN3,  'infExact','meanConst','covBranchingProcess_v5','likGauss',Data.X3,Data.Y3(:,i));         % optimise
%[ymu3 ys23 fmu3 fs23   ]= feval(@gp,hyp_pN3, 'infExact','meanConst','covBranchingProcess_v5','likGauss',Data.X3,Data.Y3(:,i),Xstar2);
 
%Shift data around
%hyp.cov  = [7;0.5;2;0.5;l1;v1;l1;v2;lg;vg]; hyp.mean = mean(Data.Y1(:,i)); hyp.lik  = 2;
%hyp_pN4 = feval(@minimize, hyp, @gp, -10000, 'infExact','meanConst','covBranchingProcess_v3','likGauss',Data.X4,Data.Y4(:,i));         % optimise
%[L4 dL4] = feval(@gp,hyp_pN4,  'infExact','meanConst','covBranchingProcess_v3','likGauss',Data.X4,Data.Y4(:,i));         % optimise
%[ymu4 ys24 fmu4 fs24   ]= feval(@gp,hyp_pN4,'infExact','meanConst','covBranchingProcess_v3','likGauss',Data.X4,Data.Y4(:,i),Xstar1);

%1->2,1->3,  Mock->DC, Mock->hrp (same as 2)
%hyp_pN5 = feval(@minimize, hyp, @gp, -10000, 'infExact','meanConst','covBranchingProcess_v4','likGauss',Data.X4,Data.Y4(:,i));         % optimise
%[L5 dL5] = feval(@gp,hyp_pN5, 'infExact','meanConst','covBranchingProcess_v4','likGauss',Data.X4,Data.Y4(:,i));         % optimise
%[ymu5 ys25 fmu5 fs25   ]= feval(@gp,hyp_pN4, 'infExact','meanConst','covBranchingProcess_v4','likGauss',Data.X4,Data.Y4(:,i),Xstar1);

%Now do joint modelling of subsets
%1,2->3,   Mock,hrp -> DC
%hyp.cov = [2;0.5;l1;v1;lg;vg];hyp.mean = mean(Data.Y1(:,i));hyp.lik = 2;
%hyp_pN6 = feval(@minimize, hyp, @gp, -10000, 'infExact','meanConst','covBranchingProcess_v5','likGauss',Data.X6,Data.Y6(:,i));         % optimise
%[L6 dL6] = feval(@gp,hyp_pN6, 'infExact','meanConst','covBranchingProcess_v5','likGauss',Data.X6,Data.Y6(:,i));         % optimise
%[ymu6 ys26 fmu6 fs26   ]= feval(@gp,hyp_pN6, 'infExact','meanConst','covBranchingProcess_v5','likGauss',Data.X6,Data.Y6(:,i),Xstar2);

% Mock-> DC,hrp
%hyp.cov = [2;0.5;l1;v1;lg;vg];hyp.mean = mean(Data.Y1(:,i));hyp.lik = 2;
%hyp_pN8 = feval(@minimize, hyp, @gp, -10000, 'infExact','meanConst','covBranchingProcess_v5','likGauss',Data.X8,Data.Y8(:,i));         % optimise
%[L8 dL8] = feval(@gp,hyp_pN8,  'infExact','meanConst','covBranchingProcess_v5','likGauss',Data.X8,Data.Y8(:,i));         % optimise
%[ymu8 ys28 fmu8 fs28   ]= feval(@gp,hyp_pN8, 'infExact','meanConst','covBranchingProcess_v5','likGauss',Data.X8,Data.Y8(:,i),Xstar3);


%Mock DC, hrp
%hyp.cov = [lg;vg]; hyp.mean = mean(Data.Y1(:,i)); hyp.lik = 2;
%hyp_pN7 = feval(@minimize, hyp, @gp, -10000, 'infExact','meanConst','covSEiso','likGauss',Data.X1(:,1),Data.Y1(:,i));         % optimise
%[L7 dL7] = feval(@gp,hyp_pN7, 'infExact','meanConst','covSEiso','likGauss',Data.X1(:,1),Data.Y1(:,i));         % optimise
%[ymu7 ys27 fmu7 fs27   ]= feval(@gp,hyp_pN7, 'infExact','meanConst','covSEiso','likGauss',Data.X1(:,1),Data.Y1(:,i),Xstar1(:,1));

%L = -[L1,L2,L3,L4,L5,L6,L7,L8];
%AIC = 2*[12,12,8,12,12,8,4,8] - 2*L;
%BIC = - 2*L + [12,12,8,12,12,8,4,8]*log(size(Data.X1,1));

%Output.L = L;
%Output.AIC = AIC;
%Output.BIC = BIC;

%Output.remove_inds = remove_inds;
%Output.ind_train = ind_train;

Ps{i,1}.remove_inds3 = remove_inds;
Ps{i,1}.ind_train3 = ind_train;


Ps{i,1}.H_2  = hyp_pN1;
Ps{i,1}.H2_2 = hyp_pN2;
Ps{i,1}.H3_2 = hyp_pN3;
Ps{i,1}.H4_2 = hyp_pN4;
%Ps{i,1}.H52 = hyp_pN5;


Ps{i,1}.L_2  = L1;
Ps{i,1}.L2_2 = L2;
Ps{i,1}.L3_2 = L3;
Ps{i,1}.L4_2 = L4;
%Ps{i,1}.L52 = L5;

Ps{i,1}.ymu1_2 = ymu1;
Ps{i,1}.ymu2_2 = ymu2;
Ps{i,1}.ymu3_2 = ymu3;
Ps{i,1}.ymu4_2 = ymu4;
%Ps{i,1}.ymu52 = ymu5;

Ps{i,1}.fmu1_2 = fmu1;
Ps{i,1}.fmu2_2 = fmu2;
Ps{i,1}.fmu3_2 = fmu3;
Ps{i,1}.fmu4_2 = fmu4;
%Ps{i,1}.fmu52 = fmu5;

Ps{i,1}.fs21_2 = fs21;
Ps{i,1}.fs22_2 = fs22;
Ps{i,1}.fs23_2 = fs23;
Ps{i,1}.fs24_2 = fs24;
%Ps{i,1}.fs252 = fs25;

Ps{i,1}.ys21_2 = ys21;
Ps{i,1}.ys22_2 = ys22;
Ps{i,1}.ys23_2 = ys23;
Ps{i,1}.ys24_2 = ys24;
%Ps{i,1}.ys252 = ys25;

%Ps{i,1} = Output;

if double(int64(i/1))==(i/1)
    save(['../results/earlydevelopment/EarlDevResults_' num2str(batchi) '_' num2str(seed) '_noprior_r.mat'],'Ps')
end

end

save(['../results/earlydevelopment/EarlDevResults_' num2str(batchi) '_' num2str(seed) '_noprior_r.mat'],'Ps')


Fin = 1;
