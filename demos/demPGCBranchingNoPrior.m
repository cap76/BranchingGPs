function [Fin] = demPGCBranchingNoPrior(batchi);

%addpath(genpath('../Pseudotime_Algorithms/Functions/gpml-matlab-v3.6-2015-07-07/'))
%addpath(('./gpml/cov'))
%addpath(genpath('./netlab3_3/'))
addpath(genpath('../'))


%load('PGClabelleddata.mat');
D1 = importdata('./AllProcesses.csv')
%D1 = importdata('/Users/christopher_penfold/Desktop/Pathway Modelling/Datasets/AllDE_P_gtp1_FC_gt10.csv')

Type = D1.data(1,:);
Sex = D1.data(2,:);
Time = D1.data(3,:);
ind = find(Time>-1);
Data.Y = log2(D1.data(4:end,ind)+1)';
Time = Time(ind);
Sex = Sex(ind);
Type = Type(ind);

z = zeros(1,length(Sex));
z(find(Type==2 & Sex==1))=ones;
z(find(Type==2 & Sex==2))=2*ones;
z(find((Type==1 | Type==0) & Sex==1))=3*ones;
z(find((Type==1 | Type==0) & Sex==2))=4*ones;
z(find(z==0))=randi([1 4],1,length(find(z==0)));

z1 = z;
z1(find(z1==2))=1;
z1(find(z1==4))=3; %All same
z2 = z;
z2(find(z2==2))=1; %PGCs same
z3 = z;
z3(find(z1==4))=3; %PGCs same

Data.X = [Time',z'];
Data.X1 = [Time',z1'];
Data.X2 = [Time',z2'];
Data.X3 = [Time',z3'];
%for i = 1:size(Data.Y,2)
%    Data.Y(:,i) = (Data.Y(:,i)-mean(Data.Y(:,i)))
%end
%keyboard

%D1 = importdata('./Pa13-Combo2.txt');

try %Try to resume analysis 
    load(['/Users/christopher_penfold/Desktop/OldBranching_GPs/PGC/PGC_new_new/PGCResults_' num2str(batchi) '_noprior_modelcomp.mat'])
    startind = length(Ps)+1;
    endind   = (batchi)*100;
catch
    startind = (batchi-1)*100 + 1;
    endind   = (batchi)*100;
end

endind = min(endind,size(Data.Y,2));

%X1 = [repmat(1:1:13,1,12); ones(1,52),2*ones(1,52),3*ones(1,52)]';
%X2 = [repmat(1:1:13,1,12); ones(1,52),ones(1,52),2*ones(1,52)]';
%X3 = [repmat(1:1:13,1,12); ones(1,52),2*ones(1,52),2*ones(1,52)]';

Xstar1 = [repmat(linspace(0,13,1000),1,4);ones(1,1000),2*ones(1,1000),3*ones(1,1000),4*ones(1,1000)]';
%Xstar2 = [repmat(linspace(0,13,50),1,3);ones(1,50),1*ones(1,50),2*ones(1,50)]';
%Xstar3 = [repmat(linspace(0,13,50),1,3);ones(1,50),2*ones(1,50),2*ones(1,50)]';

%keyboard

for i = startind:endind
    i
%Y1 = D1.data(i,:)'; %Mock hrp DC
%Y2 = D1.data(i,[1:52,2*52+1:3*52,52+1:2*52])';
l1 = log(3); l2 = log(3); lg = log(3); v1 = log(3); v2 = log(3); vg = log(std(Data.Y(:,i)));

%1->2->3. Mock->hrp->DC
hyp.cov  = [7;1.5;7;1.5;11;1.5;11;1.5;11;1.5;11;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5]; hyp.mean = mean(Data.Y(:,i)); hyp.lik  = 2;
hyp_pN1 = feval(@minimize, hyp, @gp, -10000, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X,Data.Y(:,i));         % optimise
[L1 dL1] = feval(@gp,hyp_pN1, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X,Data.Y(:,i));         % optimise
[ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X,Data.Y(:,i),Xstar1);



%1->2->3. Mock->hrp->DC
hyp.cov  = [7;1.5;7;1.5;11;1.5;11;1.5;11;1.5;11;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5]; hyp.mean = mean(Data.Y(:,i)); hyp.lik  = 2;

hyp_pN2 = feval(@minimize, hyp, @gp, -10000, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X1,Data.Y(:,i));         % optimise
[L2 dL2] = feval(@gp,hyp_pN1, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X1,Data.Y(:,i));         % optimise
[ymu2 ys22 fmu2 fs22   ]= feval(@gp,hyp_pN2, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X1,Data.Y(:,i),Xstar1);


%1->2->3. Mock->hrp->DC
hyp.cov  = [7;1.5;7;1.5;11;1.5;11;1.5;11;1.5;11;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5]; hyp.mean = mean(Data.Y(:,i)); hyp.lik  = 2;
hyp_pN3 = feval(@minimize, hyp, @gp, -10000, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X2,Data.Y(:,i));         % optimise
[L3 dL3] = feval(@gp,hyp_pN1, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X2,Data.Y(:,i));         % optimise
[ymu3 ys23 fmu3 fs23   ]= feval(@gp,hyp_pN3, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X2,Data.Y(:,i),Xstar1);


%1->2->3. Mock->hrp->DC
hyp.cov  = [7;1.5;7;1.5;11;1.5;11;1.5;11;1.5;11;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5]; hyp.mean = mean(Data.Y(:,i)); hyp.lik  = 2;
hyp_pN4 = feval(@minimize, hyp, @gp, -10000, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X3,Data.Y(:,i));         % optimise
[L4 dL4] = feval(@gp,hyp_pN4, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X3,Data.Y(:,i));         % optimise
[ymu4 ys24 fmu4 fs24   ]= feval(@gp,hyp_pN4, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X3,Data.Y(:,i),Xstar1);

k1={'covMaterniso',3};
hyp.cov  = [1.5;1.5]; hyp.mean = mean(Data.Y(:,i)); hyp.lik  = 2;
hyp_pN5 = feval(@minimize, hyp, @gp, -10000, 'infExact', 'meanConst',k1,'likGauss',Data.X3(:,1),Data.Y(:,i));         % optimise
[L5 dL5] = feval(@gp,hyp_pN5, 'infExact', 'meanConst',k1,'likGauss',Data.X3(:,1),Data.Y(:,i));         % optimise
[ymu5 ys25 fmu5 fs25   ]= feval(@gp,hyp_pN5, 'infExact', 'meanConst',k1,'likGauss',Data.X3(:,1),Data.Y(:,i),Xstar1(:,1));


%keyboard

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
Output.H = hyp_pN1;
Output.H2 = hyp_pN2;
Output.H3 = hyp_pN3;
Output.H4 = hyp_pN4;
Output.H5 = hyp_pN5;
Output.L1 = L1;
Output.L2 = L2;
Output.L3 = L3;
Output.L4 = L4;
Output.L5 = L5;
%ymu1 ys21 fmu1 fs21
%Output.H2 = hyp_pN2;
%Output.H3 = hyp_pN3;
%Output.H4 = hyp_pN4;
%Output.H5 = hyp_pN5;
%Output.H6 = hyp_pN6;
%Output.H7 = hyp_pN7;
%Output.H8 = hyp_pN8;
Output.fmu1 = fmu1;
Output.fmu2 = fmu2;
Output.fmu3 = fmu3;
Output.fmu4 = fmu4;
Output.fmu5 = fmu5;
Output.ymu1 = ymu1;
Output.fmu2 = fmu2;
Output.fmu3 = fmu3;
Output.fmu4 = fmu4;
Output.fmu5 = fmu5;
%Output.fmu6 = fmu6;
%Output.fmu7 = fmu7;
%Output.fmu8 = fmu8;
Output.fs21 = fs21;
Output.ys21 = ys21;
Output.fs22 = fs22;
Output.fs23 = fs23;
Output.fs24 = fs24;
Output.fs25 = fs25;
%Output.fs26 = fs26;
%Output.fs27 = fs27;
%Output.fs28 = fs28;

Ps{i,1} = Output;

%keyboard
if double(int64(i/1))==(i/1)
    save(['/Users/christopher_penfold/Desktop/OldBranching_GPs/PGC/PGC_new_new/PGCResults_' num2str(batchi) '_noprior_modelcomp_r.mat'],'Ps')
end

end

save(['/Users/christopher_penfold/Desktop/OldBranching_GPs/PGC/PGC_new_new/PGCResults_' num2str(batchi) '_noprior_modelcomp_r.mat'],'Ps')


Fin = 1;
