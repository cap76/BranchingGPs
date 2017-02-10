%function Output = ESCvsED(branchi);
addpath(genpath('..'))
D1 = importdata('./MarkerGenes_Naokos_Paper.xlsx');

%D1 = importdata('MarkerGenes_Naokos_Paper.xlsx');
 T = [];
 for i = 1:13
     load(['../results/primordial_germ_cells/Pseudotime/Marker_Pseudotime_1_' num2str(i) '_2E.mat'])
     T = [T,Output.Param.Store{end}.update.x(:,1)];
 end
 
 %Now pick on of the runs.
 load(['../results/primordial_germ_cells/Pseudotime/Marker_Pseudotime_1_' num2str(1) '_2E.mat'])
 ut = unique(Output.Param.Store{end}.update.origt);
 meanT = mean(T,2);
 stdT = std(T')';
 trueType = D1.data(1,:)';
 truetime = D1.data(3,:)';
 
 Type = Output.Param.Store{end}.orig.Type;
 Sex  = Output.Param.Store{end}.orig.Sex;
 infertime = Output.Param.Store{end}.update.x(:,1);
 infertype = Output.Param.Store{end}.update.x(:,2);
 
 ind1 = find( trueType==0 | (trueType==1 & Sex==1) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Female PGC
 ind2 = find( trueType==0 | (trueType==1 & Sex==2) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Male PGC
 ind3 = find( trueType==0 | (trueType==2 & Sex==1) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Female soma
 ind4 = find( trueType==0 | (trueType==2 & Sex==2) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Male soma
 
 ind5 = find((trueType==-1 | trueType==0) & Output.Param.Store{end}.update.origt~=ut(6) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2);
 ind5a = find(Output.Param.Store{end}.update.origt==ut(6) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2);
 
 ind6 = find(trueType==0 & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Early dev
 
 %Get the "ESC" time series
 ind7a = find( trueType==-1 & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %ESC
 %Range: 0.332 - 0.5352
 ind7b = find(trueType==0 & truetime~=6 & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Early dev not overlapping with ESC
 ind7 = [ind7a;ind7b];
 
% for i = 1:44
% plot(Output.Param.Store{end}.update.x(ind6,1),Output.Param.Store{end}.update.y(ind6,i),'bo')
% hold on
% plot(Output.Param.Store{end}.update.x(ind7,1)+0.01,Output.Param.Store{end}.update.y(ind7,i),'rs')
% pause
% clf
% end

%startind = (batchi-1)*1000+1;
%endind = (batchi)*1000;
%endind = min(endind,18034-3);
load(['ESCvED_1_2E.mat']) 


%D2 = importdata('GSE36552_and_GSE63818_all.csv');
D2 = importdata('./AllProcesses.csv');

im = {@infExact};                % inference method
Xstar = [linspace(0,0.5,100)',ones(100,1); linspace(0,0.5,100)',2*ones(100,1)];
for i = 2246:size(D2.data,1)-3

    x1 = Output.Param.Store{end}.update.x(ind6,1);
    x2 = Output.Param.Store{end}.update.x(ind7,1);
    y1 = D2.data(i+3,ind6)';
    y2 = D2.data(i+3,ind7)';
    
    X = [x1,ones(size(x1,1),1); x2,2*ones(size(x2,1),1)];
    Y = [y1;y2];
    
     %Joint GP (not DE). 
    k1={'covMaterniso',3}; 
    hyp.cov = [log(3);log(3)]; hyp.mean = mean(Y(:,1)); hyp.lik = log(2);
    prior.mean = {[]};  prior.cov  = {[];[]}; prior.lik = {[]};
    par1a = {'meanConst',k1,'likGauss',X(:,1),Y};
    par1b = {'meanConst',k1,'likGauss',X(:,1),Y,Xstar(:,1)};
    hyp_pN2 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L2 dL2] = feval(@gp,hyp_pN2, im, par1a{:});         % optimise
    [ymu2 ys22 fmu2 fs22   ]= feval(@gp,hyp_pN2, im, par1b{:});
    
    hyp.cov = [0.4; 2;hyp_pN2.cov;hyp_pN2.cov];
    hyp.mean = hyp_pN2.mean; hyp.lik = hyp_pN2.lik;

    par1a = {'meanConst','covBranchingProcess_2E','likGauss',X,Y};
    par1b = {'meanConst','covBranchingProcess_2E','likGauss',X,Y,Xstar};
    hyp_pN1 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L1 dL1] = feval(@gp,hyp_pN1, im, par1a{:});         % optimise
    [ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, im, par1b{:});

    %Store likelihoods
    L   = -[L1,L2];
    AIC = 2*[8,4] - 2*L;
    BIC = - 2*L + [8,4]*log(size(X,1));

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
    ESCvED{i}.hyp1 = hyp_pN1;
    ESCvED{i}.hyp2 = hyp_pN2;    

    save(['ESCvED_1_2E.mat'],'ESCvED')    
    disp(['Step ' num2str(i)])
end







 trueType = D1.data(1,:)';
 truetime = D1.data(3,:)';
 
 Type = Output.Param.Store{end}.orig.Type;
 Sex  = Output.Param.Store{end}.orig.Sex;
 infertime = Output.Param.Store{end}.update.x(:,1);
 infertype = Output.Param.Store{end}.update.x(:,2);
 
 ind1 = find( trueType==0 | (trueType==1 & Sex==1) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Female PGC
 ind2 = find( trueType==0 | (trueType==1 & Sex==2) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Male PGC
 ind3 = find( trueType==0 | (trueType==2 & Sex==1) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Female soma
 ind4 = find( trueType==0 | (trueType==2 & Sex==2) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Male soma
 
 ind5 = find((trueType==-1 | trueType==0) & Output.Param.Store{end}.update.origt~=ut(6) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2);
 ind5a = find(Output.Param.Store{end}.update.origt==ut(6) & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2);
 
 ind6 = find(trueType==0 & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2 & infertype==1); %Early dev

 s1 = randperm(length(ind2));
  ind6 = [ind6;ind2(s1(1:double(int64(length(s1/2)))))];
 
 %Get the "ESC" time series
 ind7a = find( trueType==-1 & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2 & infertype==1); %ESC
 %Range: 0.332 - 0.5352
 ind7b = find(trueType==0 & truetime~=6 & Output.Param.Store{end}.update.x(:,1)<0.9 & stdT<0.2); %Early dev not overlapping with ESC
 ind7 = [ind7a;ind7b;ind2(s1(double(int64(length(s1/2)))+1:end))];
 

im = {@infExact};                % inference method
Xstar = [linspace(0,0.5,100)',ones(100,1); linspace(0,0.5,100)',2*ones(100,1)];
for i = 2246:size(D2.data,1)-3

    x1 = Output.Param.Store{end}.update.x(ind6,1);
    x2 = Output.Param.Store{end}.update.x(ind7,1);
    y1 = D2.data(i+3,ind6)';
    y2 = D2.data(i+3,ind7)';
    
    X = [x1,ones(size(x1,1),1); x2,2*ones(size(x2,1),1)];
    Y = [y1;y2];
    
       
    %Joint GP (not DE). 
    k1={'covMaterniso',3}; 
    hyp.cov = [log(3);log(3)]; hyp.mean = mean(Y(:,1)); hyp.lik = log(2);
    prior.mean = {[]};  prior.cov  = {[];[]}; prior.lik = {[]};
    par1a = {'meanConst',k1,'likGauss',X(:,1),Y};
    par1b = {'meanConst',k1,'likGauss',X(:,1),Y,Xstar(:,1)};
    hyp_pN2 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L2 dL2] = feval(@gp,hyp_pN2, im, par1a{:});         % optimise
    [ymu2 ys22 fmu2 fs22   ]= feval(@gp,hyp_pN2, im, par1b{:});
    
    pcp1     = {@priorSmoothBox1,0,2,10};    % Gaussian prior
    pcp2     = {@priorGauss,0,1}; 
    
    hyp.cov  = [0.5;1;0.3;1;hyp_pN2.cov;hyp_pN2.cov]; hyp.mean = mean(Y(:,1)); hyp.lik  = 2;
    prior.mean = {[]};  prior.cov  = {pcp1;pcp2;pcp1;pcp2;[];[];[];[]}; prior.lik = {[]}; %prior.cov  = {pcp1p2;[];[];pcp1;[];[];[];[];[];[];[];[]}; prior.lik = {[]};
    im = {@infPrior,@infExact,prior};                % inference method
    par1a = {'meanConst','covBranchingRecombinationProcess_2C','likGauss',X,Y};
    par1b = {'meanConst','covBranchingRecombinationProcess_2C','likGauss',X,Y,Xstar};
    hyp_pN1 = feval(@minimize, hyp, @gp, -40000, im, par1a{:});         % optimise
    [L1 dL1] = gp(hyp_pN1, im, par1a{:});         % optimise
    [ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, im, par1b{:});

    
    %hyp.cov = [0.4; 2;hyp_pN2.cov;hyp_pN2.cov];
    %hyp.mean = hyp_pN2.mean; hyp.lik = hyp_pN2.lik;

    %par1a = {'meanConst','covBranchingProcess_2E','likGauss',X,Y};
    %par1b = {'meanConst','covBranchingProcess_2E','likGauss',X,Y,Xstar};
    %hyp_pN1 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    %[L1 dL1] = feval(@gp,hyp_pN1, im, par1a{:});         % optimise
    %[ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, im, par1b{:});

    %Store likelihoods
    L   = -[L1,L2];
    AIC = 2*[8,4] - 2*L;
    BIC = - 2*L + [8,4]*log(size(X,1));

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
    ESCvED{i}.hyp1 = hyp_pN1;
    ESCvED{i}.hyp2 = hyp_pN2;    

    save(['ESCvED_1_2E_recomb.mat'],'ESCvED')    
    disp(['Step ' num2str(i)])
end
