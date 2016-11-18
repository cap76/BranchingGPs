function [Fin] = demPGCBranchingNoPrior(batchi);

%Get branching process over PGC datasets. A four component system.

addpath(genpath('../'))

%addpath(genpath('../Pseudotime_Algorithms/Functions/gpml-matlab-v3.6-2015-07-07/'))
%addpath(('./gpml/cov'))
%addpath(genpath('./netlab3_3/'))
 
%Load PGC data and get into correct format
D1 = importdata('./AllProcesses.csv')
Type = D1.data(1,:);
Sex  = D1.data(2,:);
Time = D1.data(3,:);
ind  = find(Time>-1); %Ignore ESC for this test
Data.Y = log2(D1.data(4:end,ind)+1)';
Time = Time(ind);
Sex  = Sex(ind);
Type = Type(ind);

%Assign label according to group (sex and type of cell)
z = zeros(1,length(Sex));
z(find(Type==2 & Sex==1))=ones;
z(find(Type==2 & Sex==2))=2*ones;
z(find((Type==1 | Type==0) & Sex==1))=3*ones;
z(find((Type==1 | Type==0) & Sex==2))=4*ones;

%Randomly assign the pre-implantation cells a label.
z(find(z==0))=randi([1 4],1,length(find(z==0)));

Data.X = [Time',z'];

try %Try to resume analysis 
    load(['./results/primordial_germ_cells/PGCResults_' num2str(batchi) '_noprior.mat'])
    startind = length(Ps)+1;
    endind   = (batchi)*300;
catch
    startind = (batchi-1)*300 + 1;
    endind   = (batchi)*300;
end
endind = min(endind,size(Data.Y,2));

Xstar1 = [repmat(linspace(0,13,1000),1,4);ones(1,1000),2*ones(1,1000),3*ones(1,1000),4*ones(1,1000)]';

for i = startind:endind

    
    keyboard
    
%l1 = log(3); l2 = log(3); lg = log(3); v1 = log(3); v2 = log(3); vg = log(std(Data.Y(:,i)));

hyp.cov  = [7;0.5;7;0.5;11;0.5;11;0.5;11;0.5;11;0.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5;1.5]; hyp.mean = mean(Data.Y(:,i)); hyp.lik  = 2;

%hyp_pN  = feval(@minimize, hyp, @gp, -50000, @infExact, 'meanConst','covBranchingProcess_4B','likGauss',X,Y);
hyp_pN1 = feval(@minimize, hyp, @gp, -10000, @infExact, 'meanConst','covBranchingProcess_4B','likGauss',Data.X,Data.Y(:,i));         % optimise
[L1 dL1] = feval(@gp,hyp_pN1, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X,Data.Y(:,i));         % optimise
[ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, 'infExact', 'meanConst','covBranchingProcess_4B','likGauss',Data.X,Data.Y(:,i),Xstar1);

Output.H = hyp_pN1;

Output.fmu1 = fmu1;
Output.ymu1 = ymu1;
Output.fs21 = fs21;
Output.ys21 = ys21;

Ps{i,1} = Output;

if double(int64(i/10))==(i/10)
    save(['./results/primordial_germ_cells/PGCResults_' num2str(batchi) '_noprior.mat'],'Ps')
end

end

save(['./results/primordial_germ_cells/PGCResults_' num2str(batchi) '_noprior.mat'],'Ps')


Fin = 1;
