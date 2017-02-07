function Output = ESCvsED(branchi);
%addpath(genpath('../../../'))
%D1 = importdata('../../../demos/MarkerGenes_Naokos_Paper.xlsx');

%D1 = importdata('MarkerGenes_Naokos_Paper.xlsx');
 T = [];
 for i = 1:13
     load(['Marker_Pseudotime_1_' num2str(i) '_2E.mat'])
 end
 meanT = mean(T,2);
 stdT = std(T')';
 trueType = D1.data(1,:)';
 truetime = D1.data(3,:)';
 
 Type = Output.Param.Store{end}.orig.Type;
 Sex  = Output.Param.Store{end}.orig.Sex;
 infertime = Output.Param.Store{end}.update.x(:,1);
 infertype = Output.Param.Store{end}.update.x(:,2);
 
 ind1 = find( trueType==0 | (trueType==1 & Sex==1) ); %Female PGC
 ind2 = find( trueType==0 | (trueType==1 & Sex==2) ); %Male PGC
 ind3 = find( trueType==0 | (trueType==2 & Sex==1) ); %Female soma
 ind4 = find( trueType==0 | (trueType==2 & Sex==2) ); %Male soma
 
 ind5 = find((trueType==-1 | trueType==0) & Output.Param.Store{end}.update.origt~=ut(6));
 ind5a = find(Output.Param.Store{end}.update.origt==ut(6));
 
 ind6 = find(trueType==0 ); %Early dev
 
 %Get the "ESC" time series
 ind7a = find( trueType==-1); %ESC
 %Range: 0.332 - 0.5352
 ind7b = find(trueType==0 & (Output.Param.Store{end}.update.x(:,1)<0.388 | Output.Param.Store{end}.update.x(:,1)>0.5032)); %Early dev not overlapping with ESC
 ind7 = [ind7a;ind7b];
 
% for i = 1:44
% plot(Output.Param.Store{end}.update.x(ind6,1),Output.Param.Store{end}.update.y(ind6,i),'bo')
% hold on
% plot(Output.Param.Store{end}.update.x(ind7,1)+0.01,Output.Param.Store{end}.update.y(ind7,i),'rs')
% pause
% clf
% end

startind = (batchi-1)*1000+1;
endind = (batchi)*1000;
endind = min(endind,18034-3);

D2 = importdata('GSE36552_and_GSE63818_all.csv');


im = {@infExact};                % inference method
Xstar = [linspace(0,0.5,100)',ones(100,1); linspace(0,0.5,100)',2*ones(100,1)];
for i = startind:endind-3

    x1 = Output.Param.Store{end}.update.x(ind6,1);
    x2 = Output.Param.Store{end}.update.x(ind7,1);
    y1 = D2.data(i+3,ind6)';
    y2 = D2.data(i+3,ind7)';
    
    X = [x1,ones(size(x1,1),1); x2,2*ones(size(x2,1),1)];
    Y = [y1;y2];
    
    hyp.cov = [0.4; 1;log(3);log(3);log(3);log(3)];
    hyp.mean = mean(Y(:,1)); hyp.lik = log(2);

    par1a = {'meanConst','covBranchingProcess_2A','likGauss',X,Y};
    par1b = {'meanConst','covBranchingProcess_2A','likGauss',X,Y,Xstar};
    hyp_pN1 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L1 dL1] = feval(@gp,hyp_pN1, im, par1a{:});         % optimise
    [ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, im, par1b{:});

    %Joint GP (not DE). 
    hyp.cov = [log(3);log(3)]; hyp.mean = mean(Y(:,1)); hyp.lik = log(2);
    prior.mean = {[]};  prior.cov  = {[];[]}; prior.lik = {[]};
    par1a = {'meanConst','covSEiso','likGauss',X(:,1),Y};
    par1b = {'meanConst','covSEiso','likGauss',X(:,1),Y,Xstar(:,1)};
    hyp_pN2 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L2 dL2] = feval(@gp,hyp_pN2, im, par1a{:});         % optimise
    [ymu2 ys22 fmu2 fs22   ]= feval(@gp,hyp_pN2, im, par1b{:});

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

    save(['ESCvED_' num2str(batchi) '.mat'],'ESCvED')    
end
