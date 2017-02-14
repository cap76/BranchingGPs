function Output = ESCvsED(batchi);

warning off all

addpath(genpath('..'))

D1 = importdata('./MarkerGenes_Naokos_Paper.xlsx');

T = [];
for i = 1:13
    load(['../results/primordial_germ_cells/Pseudotime/Marker_Pseudotime_1_' num2str(i) '_2E.mat'])
    T = [T,Output.Param.Store{end}.update.x(:,1)];
end
 
%Now pick one of the runs.
load(['../results/primordial_germ_cells/Pseudotime/Marker_Pseudotime_1_' num2str(6) '_2E.mat'])
ut = unique(Output.Param.Store{end}.update.origt);
meanT = mean(T,2);
stdT = std(T')';
trueType = D1.data(1,:)';
truetime = D1.data(3,:)';
 
trueType = D1.data(1,:)';
truetime = D1.data(3,:)';
 
Type = Output.Param.Store{end}.orig.Type;
Sex  = Output.Param.Store{end}.orig.Sex;
infertime = Output.Param.Store{end}.update.x(:,1);
infertype = Output.Param.Store{end}.update.x(:,2);
 
ind1 = find( (trueType==1 & Sex==1) & Output.Param.Store{end}.update.x(:,1)<0.9 ); %Female PGC
ind2 = find( (trueType==1 & Sex==2) & Output.Param.Store{end}.update.x(:,1)<0.9 ); %Male PGC
ind3 = find( (trueType==2 & Sex==1) & Output.Param.Store{end}.update.x(:,1)<0.9 ); %Female soma
ind4 = find( (trueType==2 & Sex==2) & Output.Param.Store{end}.update.x(:,1)<0.9 ); %Male soma
 
ind5 = find((trueType==-1 | trueType==0) & Output.Param.Store{end}.update.origt~=ut(6) & Output.Param.Store{end}.update.x(:,1)<0.9 );
ind5a = find(Output.Param.Store{end}.update.origt==ut(6) & Output.Param.Store{end}.update.x(:,1)<0.9 );
 
%Ge the early development dataset
ind6 = find(trueType==0 & Output.Param.Store{end}.update.x(:,1)<0.9 ); %Early dev
s1 = randperm(length(ind2));
%Add some PGC data 
ind6 = [ind6;ind2(s1(1:double(int64(length(s1)/2))))];
 
%Get the "ESC" time series
ind7a = find(trueType==-1 & Output.Param.Store{end}.update.x(:,1)<0.9); %ESC
%Get pre-blastocyst
ind7b = find(trueType==0 & truetime~=6 & Output.Param.Store{end}.update.x(:,1)<0.9 ); %Early dev not overlapping with ESC
%Add some PGC
ind7 = [ind7a;ind7b;ind2(s1(double(int64(length(s1)/2))+1:end))];

Xstar = [linspace(0,0.8,100)',ones(100,1); linspace(0,0.8,100)',2*ones(100,1)];

%D2 = importdata('AllProcesses.csv')
D2 = D1;
%D2 = importdata('GSE36552_and_GSE63818_all.csv');
 
genes = D2.textdata(4:end,1);

startind = (batchi-1)*1000+1;
endind = (batchi)*1000;
endind = min(endind,size(D2.data,1)-3);

genes = D2.textdata(6:end,1);

for i = 1:size(D2.data,1)-3

    try
        
    x1 = Output.Param.Store{end}.update.x(ind6,1);
    x2 = Output.Param.Store{end}.update.x(ind7,1);
    y1 = log2(D2.data(i+3,ind6)'+1);
    y2 = log2(D2.data(i+3,ind7)'+1);
    
    X = [x1,ones(size(x1,1),1); x2,2*ones(size(x2,1),1)];
    Y = [y1;y2];
           
    %Joint GP (not DE). 
    pLS     = {@priorGauss,0.5,1};
    pN      = {@priorClamped}; 
    
    hyp.cov = [log(3);log(3)]; hyp.mean = mean(Y(:,1)); hyp.lik = log(2);
    prior.mean = {[]};  prior.cov  = {pLS;[]}; prior.lik = {[]};
    im = {@infPrior,@infExact,prior}; 
    par1a = {'meanConst',@covSEiso,'likGauss',X(:,1),Y};
    par1b = {'meanConst',@covSEiso,'likGauss',X(:,1),Y,Xstar(:,1)};
    hyp_pN2 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L2 dL2] = feval(@gp,hyp_pN2, im, par1a{:});         % optimise
    [ymu2 ys22 fmu2 fs22   ]= feval(@gp,hyp_pN2, im, par1b{:});
    
    pcp1     = {@priorSmoothBox1,0,2,10};    % Gaussian prior
    pcp2     = {@priorGauss,2,2}; 
    
    hyp.cov = [0.4; 3; hyp_pN2.cov;hyp_pN2.cov];    
    hyp.mean = hyp_pN2.mean; hyp.lik = hyp_pN2.lik;
    prior.mean = {[]};  prior.cov  = {pcp1;pcp2;pLS;[];pLS;[]}; prior.lik = {pN};
    im = {@infPrior,@infExact,prior};     
    par1a = {'meanConst','covBranchingProcess_2A','likGauss',X,Y};
    par1b = {'meanConst','covBranchingProcess_2A','likGauss',X,Y,Xstar};
    hyp_pN3 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L3 dL3] = feval(@gp,hyp_pN3, im, par1a{:});         % optimise
    [ymu3 ys23 fmu3 fs23   ]= feval(@gp,hyp_pN3, im, par1b{:});
    
    hyp.cov  = [hyp_pN3.cov(1)+0.2;hyp_pN3.cov(2);hyp_pN3.cov(1:2);hyp_pN3.cov(3:end)]; hyp.mean = mean(Y(:,1)); hyp.lik  = log(2);
    prior.mean = {[]};  prior.cov  = {pcp1;pcp2;pcp1;pcp2;pLS;[];pLS;[]}; prior.lik = {pN}; %prior.cov  = {pcp1p2;[];[];pcp1;[];[];[];[];[];[];[];[]}; prior.lik = {[]};
    im = {@infPrior,@infExact,prior};                % inference method
    par1a = {'meanConst','covBranchingRecombinationProcess_2C','likGauss',X,Y};
    par1b = {'meanConst','covBranchingRecombinationProcess_2C','likGauss',X,Y,Xstar};
    hyp_pN1 = feval(@minimize, hyp, @gp, -40000, im, par1a{:});         % optimise
    [L1 dL1] = gp(hyp_pN1, im, par1a{:});         % optimise
    [ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, im, par1b{:});
       
    %Store likelihoods
    L   = -[L1,L2,L3];
    AIC = 2*[8,4,6] - 2*L;
    BIC = - 2*L + [8,4,6]*log(size(X,1));

    ESCvED{i}.gene = genes{i};
    
    ESCvED{i}.L = L;
    ESCvED{i}.AIC = AIC;
    ESCvED{i}.BIC = BIC;
    ESCvED{i}.fmu1 = fmu1;
    ESCvED{i}.fs21 = fs21;    
    ESCvED{i}.ymu1 = ymu1;
    ESCvED{i}.ys21 = ys21;
    
    ESCvED{i}.fmu2 = fmu2;
    ESCvED{i}.fs22 = fs22;    
    ESCvED{i}.ymu2 = ymu2;
    ESCvED{i}.ys22 = ys22;    

    ESCvED{i}.fmu3 = fmu3;
    ESCvED{i}.fs23 = fs23;    
    ESCvED{i}.ymu3 = ymu3;
    ESCvED{i}.ys23 = ys23;    

    
    ESCvED{i}.hyp1 = hyp_pN1;
    ESCvED{i}.hyp2 = hyp_pN2;    
    ESCvED{i}.hyp3 = hyp_pN3;    
    

    save(['ESCvED_1_2E_' num2str(batchi) '_NaokoMarkers.mat'],'ESCvED')    
    disp(['Step ' num2str(i)])
    catch
    end
end

Output = 1;