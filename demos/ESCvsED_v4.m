function Output = ESCvsED_v4(batchi);

warning off all

addpath(genpath('..'))

D1 = importdata('./MarkerGenes_Naokos_Paper.xlsx');

T = [];
for i = 1:3
    load(['../results/primordial_germ_cells/Pseudotime/Marker_Pseudotime_' num2str(i) '_' num2str(6) '_2Ec.mat'])
    T = [T,Output.Param.Store{end}.update.x(:,1)];
end
 
%Now pick one of the runs.
load(['../results/primordial_germ_cells/Pseudotime/Marker_Pseudotime_1_6_2Ec.mat'])
ut    = unique(Output.Param.Store{end}.update.origt);
meanT = mean(T,2);
stdT  = std(T')';
trueType = D1.data(1,:)';
truetime = D1.data(3,:)';
 
trueType = D1.data(1,:)';
truetime = D1.data(3,:)';
 
Type = Output.Param.Store{end}.orig.Type;
Sex  = Output.Param.Store{end}.orig.Sex;
infertime = Output.Param.Store{end}.update.x(:,1);
infertype = Output.Param.Store{end}.update.x(:,2);
 
ind1 = find( (trueType==1 & Sex==1) ); %Female PGC

ind2 = find( (trueType==1 & Sex==2) ); %Male PGC
    r2 = randperm(length(ind2));
    ind2a = ind2(r2(1:double(int64(length(r2)/2))));
    ind2b = ind2(r2(double(int64(length(r2)/2))+1:end));    
    
ind3 = find( (trueType==2) ); %Soma

ind5 = find( truetime<6   ); %Pre-blastocyst
       r5 = randperm(length(ind5));
    ind5a = ind5(r5(1:double(int64(length(r5)/4))));
    ind5b = ind5(r5(double(int64(length(r5)/4))+1:2*double(int64(length(r5)/4))));    
    ind5c = ind5(r5(2*double(int64(length(r5)/4))+1:3*double(int64(length(r5)/4))));
    ind5d = ind5(r5(3*double(int64(length(r5)/4))+1:end));
    
ind6 = find( truetime==6  ); %Blastocyst
    r6 = randperm(length(ind6));
    ind6a = ind6(r6(1:double(int64(length(r6)/3))));
    ind6b = ind6(r6(double(int64(length(r6)/3))+1:2*double(int64(length(r6)/3))));    
    ind6c = ind6(r6(2*double(int64(length(r6)/4))+1:end));

ind7  = find( truetime==-1 ); %ESC
ind8  = [ind5a; ind6a; ind2a];
ind9  = [ind5d; ind7;  ind2b];
ind10 = [ind5c; ind6b; ind1];
ind11 = [ind5d; ind6c; ind3];

Xstar = [linspace(0,1.1,1000)',zeros(1000,1); linspace(0,1.1,1000)',ones(1000,1); linspace(0,1.1,1000)',2*ones(1000,1); linspace(0,1.1,1000)',3*ones(1000,1); linspace(0,1.1,1000)',4*ones(1000,1)];

D2 = D1;
%D2 = importdata('GSE36552_and_GSE63818_all.csv');
 
genes = D2.textdata(4:end,1);

startind = (batchi-1)*1000+1;
endind = (batchi)*1000;
endind = min(endind,size(D2.data,1));

genes = D2.textdata(6:end,1);
for i = startind:endind
    try
        
    cT1 = {@priorGamma,9,.05};   %0.4, .35
    cT2 = {@priorGamma,11,.05};  %0.55, .5
    cT3 = {@priorGamma,16,.05};  %0.8, 0.75         
    pT  = {@priorGauss,3,1}; 
    
    x1 = Output.Param.Store{end}.update.x(ind8,1);
    x2 = Output.Param.Store{end}.update.x(ind9,1);    
    x3 = Output.Param.Store{end}.update.x(ind10,1); 
    x4 = Output.Param.Store{end}.update.x(ind11,1);    
    
    y1 = log2(D2.data(i+3,ind8)'+1);
    y2 = log2(D2.data(i+3,ind9)'+1);
    y3 = log2(D2.data(i+3,ind10)'+1);
    y4 = log2(D2.data(i+3,ind11)'+1);
        
    X  = [x1,0*ones(size(x1,1),1); x2,ones(size(x2,1),1); x3,2*ones(size(x3,1),1); x4,3*ones(size(x4,1),1)];  %Everything different      
    
    Y = [y1;y2;y3;y4];
           
    %Joint GP (not DE). 
    pLS     = {@priorGauss,0.5,1}; %Prior over length scale.
    pN      = {@priorGauss,0.5,1};
    k1={'covMaterniso',3};
      
    
    %First fit base process (MPGC)
    hyp.cov = [log(3);log(3)]; hyp.mean = mean(Y(:,1)); hyp.lik = log(2);
    prior.mean = {[]};  prior.cov  = {pLS;[]}; prior.lik = {[]};
    im = {@infPrior,@infExact,prior}; 
    par1a = {'meanConst',k1,'likGauss',X(find(X(:,2)==0),1),Y(find(X(:,2)==0))};
    par1b = {'meanConst',k1,'likGauss',X(find(X(:,2)==0),1),Y(find(X(:,2)==0)),Xstar(:,1)};
    par1c = {'meanConst',k1,'likGauss',X(find(X(:,2)==0),1),Y(find(X(:,2)==0)),X(find(X(:,2)==1),1)}; %Predict at female locations    
    par1d = {'meanConst',k1,'likGauss',X(find(X(:,2)==0),1),Y(find(X(:,2)==0)),X(find(X(:,2)==2),1)}; %Predict at ESC locations        
    par1e = {'meanConst',k1,'likGauss',X(find(X(:,2)==0),1),Y(find(X(:,2)==0)),X(find(X(:,2)==3),1)}; %Predict at ESC locations            
    hyp_pN1 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L1 dL1] = feval(@gp,hyp_pN1, im, par1a{:});         % optimise
    
    %
    [ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, im, par1b{:});
    [ymu1c ys21c fmu1c fs21c   ]= feval(@gp,hyp_pN1, im, par1c{:});
    [ymu1d ys21d fmu1d fs21d   ]= feval(@gp,hyp_pN1, im, par1d{:});
    [ymu1e ys21e fmu1e fs21e   ]= feval(@gp,hyp_pN1, im, par1e{:});

       
    %Now initialise hyperparameters over ESC.
    clear prior hyp
    k1={@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}};
     hyp.cov = [0.4;4;log(3);log(3)]; hyp.mean = []; hyp.lik = hyp_pN1.lik;
    %prior.mean = {[]};  
    prior.cov  = {cT2;pT;pLS;[]}; prior.lik = {pLS};
    im = {@infPrior,@infExact,prior}; 
    par1a = {'meanZero',k1,'likGauss',X(find(X(:,2)==1),1),Y(find(X(:,2)==1))-ymu1c};
    par1b = {'meanZero',k1,'likGauss',X(find(X(:,2)==1),1),Y(find(X(:,2)==1))-ymu1c,Xstar(:,1)};
    hyp_pN2 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L2 dL2] = feval(@gp,hyp_pN2, im, par1a{:});         % optimise    
    [ymu2 ys22 fmu2 fs22   ]= feval(@gp,hyp_pN2, im, par1b{:});    
    
    %Now initialise hyperparameters FPGC.
    hyp.cov = [0.4;4;log(3);log(3)]; hyp.mean = []; hyp.lik =  hyp_pN1.lik;
    %prior.mean = {[]};  
    prior.cov  = {cT3;pT;pLS;[]}; prior.lik = { pLS};
    
    im = {@infPrior,@infExact,prior}; 
    par1a = {'meanZero',k1,'likGauss',X(find(X(:,2)==2),1),Y(find(X(:,2)==2))-ymu1d};
    par1b = {'meanZero',k1,'likGauss',X(find(X(:,2)==2),1),Y(find(X(:,2)==2))-ymu1d,Xstar(:,1)};
    hyp_pN3 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L3 dL3] = feval(@gp,hyp_pN3, im, par1a{:});         % optimise    
    [ymu3 ys23 fmu3 fs23   ]= feval(@gp,hyp_pN3, im, par1b{:});    

    
    %Now initialise hyperparameters soma.
    hyp.cov = [0.4;4;log(3);log(3)]; hyp.mean = []; hyp.lik = hyp_pN1.lik;
    %prior.mean = {[]};  
    prior.cov  = {cT1;pT;pLS;[]}; prior.lik = { pLS};
    
    im = {@infPrior,@infExact,prior}; 
    par1a = {'meanZero',k1,'likGauss',X(find(X(:,2)==3),1),Y(find(X(:,2)==3))-ymu1e};
    par1b = {'meanZero',k1,'likGauss',X(find(X(:,2)==3),1),Y(find(X(:,2)==3))-ymu1e,Xstar(:,1)};
    hyp_pN4 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L4 dL4] = feval(@gp,hyp_pN4, im, par1a{:});         % optimise        
    [ymu4 ys24 fmu4 fs24   ]= feval(@gp,hyp_pN4, im, par1b{:});                
    
    %Now that we've initialised everything, let's start with more complex branching patterns.
    pC = {@priorClamped};       %Mean
    %pLS     = {@priorGauss,0.5,1};

            
    clear prior hyp
    %1) All different.    
    hyp.cov = [hyp_pN2.cov(1)+0.1;hyp_pN2.cov(2);hyp_pN2.cov(1:2);hyp_pN3.cov(1:2);hyp_pN4.cov(1:2);  hyp_pN2.cov(3:4);hyp_pN3.cov(3:4);hyp_pN4.cov(3:4);  hyp_pN1.cov(1:2)];
    hyp.mean = mean(Y(:,1)); hyp.lik = hyp_pN1.lik;    
    prior.mean = {[]};  prior.cov  = {cT2;pT;cT1;pT;cT1;pT;cT3;pT;pLS;[];pLS;[];pLS;[];pLS;[]}; prior.lik = {pN};
        
    %prior.mean = {[]};  prior.cov  = {[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[]}; prior.lik = {pLS};    
    im = {@infPrior,@infExact,prior};     
    par1a = {'meanConst','covBranchingProcess_5Rec','likGauss',X,Y};
    par1b = {'meanConst','covBranchingProcess_5Rec','likGauss',X,Y,Xstar};
    hyp_pN5 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L5 dL5] = feval(@gp,hyp_pN5, im, par1a{:});         % optimise
    [ymu5 ys25 fmu5 fs25   ]= feval(@gp,hyp_pN5, im, par1b{:});
        
    %2) Only soma is diff. 
    X1 = [x1,0*ones(size(x1,1),1); x2,0*ones(size(x2,1),1); x3,0*ones(size(x3,1),1); x4,3*ones(size(x4,1),1)]; %Only soma different    

    hyp.cov = [4;hyp_pN2.cov(2);2;hyp_pN2.cov(2); 2; hyp_pN3.cov(2); hyp_pN4.cov(1:2);  hyp_pN2.cov(3:4);hyp_pN3.cov(3:4);hyp_pN4.cov(3:4);  hyp_pN1.cov(1:2)];
    hyp.mean = mean(Y(:,1)); hyp.lik = hyp_pN1.lik;    
    prior.mean = {[]};  prior.cov  = {pC;pC;pC;pC;pC;pC;cT3;pT; pLS;[]; pC;pC; pC;pC; pLS;[]}; prior.lik = {pN};
    im = {@infPrior,@infExact,prior};     
    par1a = {'meanConst','covBranchingProcess_5Rec','likGauss',X1,Y};
    par1b = {'meanConst','covBranchingProcess_5Rec','likGauss',X1,Y,Xstar};
    hyp_pN6 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L6 dL6] = feval(@gp,hyp_pN5, im, par1a{:});         % optimise
    [ymu6 ys26 fmu6 fs26   ]= feval(@gp,hyp_pN6, im, par1b{:});    

    %3) Only FPGC is diff ... junk    
    X2 = [x1,0*ones(size(x1,1),1); x2,0*ones(size(x2,1),1); x3,2*ones(size(x3,1),1); x4,0*ones(size(x4,1),1)];       
    hyp.cov = [4;hyp_pN2.cov(2);2;hyp_pN2.cov(2); hyp_pN3.cov(1:2); 2; hyp_pN4.cov(2);  hyp_pN2.cov(3:4);hyp_pN3.cov(3:4);hyp_pN4.cov(3:4);  hyp_pN1.cov(1:2)];
    hyp.mean = mean(Y(:,1)); hyp.lik = hyp_pN1.lik;    
    prior.mean = {[]};  prior.cov  = {pC;pC;pC;pC; cT1;pT; pC;pC; pLS;[]; pC;pC; pLS;[]; pC;pC}; prior.lik = {pN};
    im = {@infPrior,@infExact,prior};     
    par1a = {'meanConst','covBranchingProcess_5Rec','likGauss',X2,Y};
    par1b = {'meanConst','covBranchingProcess_5Rec','likGauss',X2,Y,Xstar};
    hyp_pN7 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L7 dL7] = feval(@gp,hyp_pN7, im, par1a{:});         % optimise
    [ymu7 ys27 fmu7 fs27   ]= feval(@gp,hyp_pN7, im, par1b{:});        

    %4) MPG and ESC same, all others different  
    X3 = [x1,0*ones(size(x1,1),1); x2,0*ones(size(x2,1),1); x3,2*ones(size(x3,1),1); x4,3*ones(size(x4,1),1)]; %Only soma different         
    hyp.cov = [4;hyp_pN2.cov(2);2;hyp_pN2.cov(2);hyp_pN3.cov(1:2); hyp_pN4.cov(1:2);  hyp_pN2.cov(3:4);hyp_pN3.cov(3:4);hyp_pN4.cov(3:4);  hyp_pN1.cov(1:2)];
    hyp.mean = mean(Y(:,1)); hyp.lik = hyp_pN1.lik;    
    prior.mean = {[]};  prior.cov  = {pC;pC;pC;pC;cT1;pT;cT3;pT; pLS;[]; pC;pC; pLS;[]; pLS;[]}; prior.lik = {pN};
    im = {@infPrior,@infExact,prior};     
    par1a = {'meanConst','covBranchingProcess_5Rec','likGauss',X3,Y};
    par1b = {'meanConst','covBranchingProcess_5Rec','likGauss',X3,Y,Xstar};
    hyp_pN8 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L8 dL8] = feval(@gp,hyp_pN8, im, par1a{:});         % optimise
    [ymu8 ys28 fmu8 fs28   ]= feval(@gp,hyp_pN8, im, par1b{:});        
    
    %5) Soma and ESC are different   
    X4 = [x1,0*ones(size(x1,1),1); x2,1*ones(size(x2,1),1); x3,0*ones(size(x3,1),1); x4,3*ones(size(x4,1),1)]; %Only soma different         
    hyp.cov = [hyp_pN2.cov(1)+0.1;hyp_pN2.cov(2);hyp_pN2.cov(1:2); 2;hyp_pN3.cov(2); hyp_pN4.cov(1:2);  hyp_pN2.cov(3:4);hyp_pN3.cov(3:4);hyp_pN4.cov(3:4);  hyp_pN1.cov(1:2)];
    hyp.mean = mean(Y(:,1)); hyp.lik = hyp_pN1.lik;    
    prior.mean = {[]};  prior.cov  = {cT2;pT;cT1;pT;pC;pC;cT3;pT; pLS;[]; pLS;[]; pC;pC; pLS;[]}; prior.lik = {pN};
    im = {@infPrior,@infExact,prior};     
    par1a = {'meanConst','covBranchingProcess_5Rec','likGauss',X4,Y};
    par1b = {'meanConst','covBranchingProcess_5Rec','likGauss',X4,Y,Xstar};
    hyp_pN9 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L9 dL9] = feval(@gp,hyp_pN9, im, par1a{:});         % optimise
    [ymu9 ys29 fmu9 fs29   ]= feval(@gp,hyp_pN9, im, par1b{:});        


    %6) PGC and ESC are different 
    X5 = [x1,0*ones(size(x1,1),1); x2,1*ones(size(x2,1),1); x3,2*ones(size(x3,1),1); x4,0*ones(size(x4,1),1)]; %Only soma different         
    hyp.cov = [hyp_pN2.cov(1)+0.1;hyp_pN2.cov(2);hyp_pN2.cov(1:2); hyp_pN3.cov(1:2); 2;hyp_pN4.cov(2);  hyp_pN2.cov(3:4); hyp_pN3.cov(3:4);hyp_pN4.cov(3:4);  hyp_pN1.cov(1:2)];
    hyp.mean = mean(Y(:,1)); hyp.lik = hyp_pN1.lik;    
    prior.mean = {[]};  prior.cov  =  {cT2;pT;cT1;pT;cT1;pT;pC;pC; pLS;[]; pLS;[]; pLS;[];pC;pC}; prior.lik = {pN};
    im = {@infPrior,@infExact,prior};     
    par1a = {'meanConst','covBranchingProcess_5Rec','likGauss',X5,Y};
    par1b = {'meanConst','covBranchingProcess_5Rec','likGauss',X5,Y,Xstar};
    hyp_pN10 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L10 dL10] = feval(@gp,hyp_pN10, im, par1a{:});         % optimise
    [ymu10 ys210 fmu10 fs210   ]= feval(@gp,hyp_pN10, im, par1b{:});     
    


    %7) Nothing different     
    X6 = [x1,0*ones(size(x1,1),1); x2,0*ones(size(x2,1),1); x3,0*ones(size(x3,1),1); x4,0*ones(size(x4,1),1)]; %Only soma different         
    hyp.cov = [4;hyp_pN2.cov(2);2;hyp_pN2.cov(2); 2; hyp_pN3.cov(2); 2; hyp_pN4.cov(2);  hyp_pN2.cov(3:4); hyp_pN3.cov(3:4); hyp_pN4.cov(3:4);  hyp_pN1.cov(1:2)];
    hyp.mean = mean(Y(:,1)); hyp.lik = hyp_pN1.lik;    
    prior.mean = {[]};  prior.cov  =  {pC;pC;pC;pC;pC;pC;pC;pC;pLS;[];pC;pC;pC;pC;pC;pC}; prior.lik = {pN};
    im = {@infPrior,@infExact,prior};     
    par1a = {'meanConst','covBranchingProcess_5Rec','likGauss',X6,Y};
    par1b = {'meanConst','covBranchingProcess_5Rec','likGauss',X6,Y,Xstar};
    hyp_pN11 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L11 dL11] = feval(@gp,hyp_pN11, im, par1a{:});         % optimise
    [ymu11 ys211 fmu11 fs211   ]= feval(@gp,hyp_pN11, im, par1b{:});         
    

    %Store likelihoods
    L   = -[L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11];

    BIC =  -2*L(5:end) + [18,   8,8,12, 14, 14, 4]*log(length(Y));
    
    ESCvED{i}.gene = genes{i};
    ESCvED{i}.L = L;
    
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
    
    ESCvED{i}.fmu4 = fmu4;
    ESCvED{i}.fs24 = fs24;    
    ESCvED{i}.ymu4 = ymu4;
    ESCvED{i}.ys24 = ys24;        

    ESCvED{i}.fmu5 = fmu5;
    ESCvED{i}.fs25 = fs25;    
    ESCvED{i}.ymu5 = ymu5;
    ESCvED{i}.ys25 = ys25;    
    
    ESCvED{i}.fmu6 = fmu6;
    ESCvED{i}.fs26 = fs26;    
    ESCvED{i}.ymu6 = ymu6;
    ESCvED{i}.ys26 = ys26;    

    ESCvED{i}.fmu7 = fmu7;
    ESCvED{i}.fs27 = fs27;    
    ESCvED{i}.ymu7 = ymu7;
    ESCvED{i}.ys27 = ys27;        

    ESCvED{i}.fmu8 = fmu8;
    ESCvED{i}.fs28 = fs28;    
    ESCvED{i}.ymu8 = ymu8;
    ESCvED{i}.ys28 = ys28;        
    
    ESCvED{i}.fmu9 = fmu9;
    ESCvED{i}.fs29 = fs29;    
    ESCvED{i}.ymu9 = ymu9;
    ESCvED{i}.ys29 = ys29;  
    
    ESCvED{i}.fmu10 = fmu10;
    ESCvED{i}.fs210 = fs210;    
    ESCvED{i}.ymu10 = ymu10;
    ESCvED{i}.ys210 = ys210;  
    
    ESCvED{i}.fmu11 = fmu11;
    ESCvED{i}.fs211 = fs211;    
    ESCvED{i}.ymu11 = ymu11;
    ESCvED{i}.ys211 = ys211;      
    
    ESCvED{i}.hyp1 = hyp_pN1;
    ESCvED{i}.hyp2 = hyp_pN2;    
    ESCvED{i}.hyp3 = hyp_pN3;    
    ESCvED{i}.hyp4 = hyp_pN4;    
    ESCvED{i}.hyp5 = hyp_pN5;    
    ESCvED{i}.hyp6 = hyp_pN6;    
    ESCvED{i}.hyp7 = hyp_pN7;    
    ESCvED{i}.hyp8 = hyp_pN8;    
    ESCvED{i}.hyp9 = hyp_pN9;    
    ESCvED{i}.hyp10 = hyp_pN10;    
    ESCvED{i}.hyp11 = hyp_pN11;               

    save(['ESCvED_1_2E_' num2str(batchi) '_AllBranching_Matern_NaokoMarkers.mat'],'ESCvED')    
    disp(['Step ' num2str(i)])
    catch
        ESCvED{i} = [];
end

end

Output = 1;