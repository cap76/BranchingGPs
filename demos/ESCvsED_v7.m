function Output = ESCvsED_v7(batchi);

warning off all

addpath(genpath('..'))
rng(1)
D1 = importdata('GSE36552_and_GSE63818_all.csv');
%D1 = importdata('MarkerGenes_Naokos_Paper.xlsx');
%D1 = importdata('../../../demos/MarkerGenes_Naokos_Paper.xlsx');
indkeep = find(D1.data(3,:)~=6);
D1.data = D1.data(:,indkeep);
indkeep = find(D1.data(2,:)~=1);
D1.data  = D1.data(:,indkeep);


%T = [];
%for i = 1:3
%    load(['../results/primordial_germ_cells/Pseudotime/ExtendedRun/Marker_Pseudotime_' num2str(i) '_' num2str(6) '_2E_final.mat'])
%    T = [T,Output.Param.Store{299}.update.x(:,1)];
%end
%load(['~/Desktop/BranchingGPs/results/primordial_germ_cells/Pseudotime/Marker_Pseudotime_' num2str(2) '_' num2str(6) '_3Cr.mat'])                        
%load(['~/Desktop/BranchingGPs/results/primordial_germ_cells/Pseudotime/3Comp/Marker_Pseudotime_' num2str(1) '_' num2str(6.1) '_3Cr.mat'],'Output')
load(['Marker_Pseudotime_' num2str(1) '_' num2str(6.1) '_3Cr.mat'])

T = [Output.Param.Store{end}.update.x(:,1)];

%Now pick one of the runs.O
ut       = unique(Output.Param.Store{299}.update.origt);
meanT    = mean(T,2);
stdT     = std(T')';
trueType = D1.data(1,:)';
truetime = D1.data(3,:)';
  
Type      = Output.Param.Store{end}.orig.Type;
Sex       = Output.Param.Store{end}.orig.Sex;
infertime = Output.Param.Store{end}.update.x(:,1);
infertype = Output.Param.Store{end}.update.x(:,2);

Xstar = [linspace(0,1.1,1000)',ones(1000,1); linspace(0,1.1,1000)',2*ones(1000,1); linspace(0,1.1,1000)',3*ones(1000,1)];
Xstar2 = [linspace(0,1.1,1000)',ones(1000,1); linspace(0,1.1,1000)',ones(1000,1); linspace(0,1.1,1000)',2*ones(1000,1)];

D2 = D1; 
%D2 = importdata('GSE36552_and_GSE63818_all.csv');
 
%genes = D2.textdata(4:end,1);
genes = D2.textdata(6:end,1);

startind = (batchi-1)*100+1;
endind   = (batchi)*100;
endind   = min(endind,size(D2.data,1));
ind8  = find(Output.Param.Store{end}.update.x(:,2)==1);
ind9  = find(Output.Param.Store{end}.update.x(:,2)==2);
ind10 = find(Output.Param.Store{end}.update.x(:,2)==3);

%genes = D2.textdata(6:end,1);
for i = startind:endind

    try
        
    %Priors for the branch times.    
    cT1 = {@priorGamma,9,.05};   %0.4,      %Blast vs ESC
    cT2 = {@priorGamma,11,.05};  %0.55,     %PGC vs soma
    cT3 = {@priorGamma,13,.05};  %0.65,     %ESC vs PGC
    pT  = {@priorGauss,3,1}; 
    %Priors for joint GP (not DE). 
    pLS = {@priorGauss,0.5,1}; %Prior over length scale.
    pN  = {@priorGauss,0.5,1};    
    %Clamped priors
    pC  = {@priorClamped};
        
    %Observations
    x1 = Output.Param.Store{end}.update.x(ind8,1);
    x2 = Output.Param.Store{end}.update.x(ind9,1);    
    x3 = Output.Param.Store{end}.update.x(ind10,1);     
    y1 = log2(D2.data(i+3,ind8)'+1);
    y2 = log2(D2.data(i+3,ind9)'+1);
    y3 = log2(D2.data(i+3,ind10)'+1);
    %Get in correct formatt    
    X  = [x1,1*ones(size(x1,1),1); x2,2*ones(size(x2,1),1); x3,3*ones(size(x3,1),1)];  %Everything different          
    Y  = [y1;y2;y3]; 
                   
    
    %COVARIANCE FUNCTIONS
    k1 = {'covMaterniso',3};
    k2 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}};               
    k3 = {@covChangePointMultiD,{1,{@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}},@covZero}};           
    k4 = 'covBranchingProcess_3Rec3';
    k5 = 'covBranchingProcess_2Rec3';
    
    %%MODULE 1: NIAVELY INIT. HYPERPARAMETERS BASED ON PAIRWISE DATASETS
    
    %First fit base process (soma)
    hyp.cov    = [log(3);log(3)]; hyp.mean = mean(Y(:,1)); hyp.lik = log([4;4]);
    prior.mean = {[]}; prior.cov  = {pLS;[]}; prior.lik  = {[];[]};
    im    = {@infPrior,@infLaplace,prior}; 
    par1a = {'meanConst',k1,'likT',X(find(X(:,2)==1),1),Y(find(X(:,2)==1))};
    par1b = {'meanConst',k1,'likT',X(find(X(:,2)==1),1),Y(find(X(:,2)==1)),Xstar(:,1)};    
    par1c = {'meanConst',k1,'likT',X(find(X(:,2)==1),1),Y(find(X(:,2)==1)),X(find(X(:,2)==1),1)}; %Predict at soma    
    par1d = {'meanConst',k1,'likT',X(find(X(:,2)==1),1),Y(find(X(:,2)==1)),X(find(X(:,2)==2),1)}; %Predict at PGC locations        
    hyp_pN1  = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         
    [L1 dL1] = feval(@gp,hyp_pN1, im, par1a{:});
    [ymu1  ys21  fmu1  fs21    ]= feval(@gp,hyp_pN1, im, par1b{:}); %Predict @ test locations
    [ymu1c ys21c fmu1c fs21c   ]= feval(@gp,hyp_pN1, im, par1c{:}); %Predict soma @ soma locations
    [ymu1d ys21d fmu1d fs21d   ]= feval(@gp,hyp_pN1, im, par1d{:}); %Predict soma @ PGC locations
    
    %Now initialise hyperparameters over PGC.
    clear prior hyp
    hyp.cov = [0.4; 4; log(3);log(3)]; hyp.mean = []; hyp.lik = log([4;4]);
    prior.cov  = {cT2;pT;pLS;[]}; prior.lik = {[];[]};
    im = {@infPrior,@infVB,prior}; 
    par1a = {'meanZero',k2,'likT',X(find(X(:,2)==2),1),Y(find(X(:,2)==2))-ymu1d};
    par1b = {'meanZero',k2,'likT',X(find(X(:,2)==2),1),Y(find(X(:,2)==2))-ymu1d,Xstar(:,1)};    
    hyp_pN2  = feval(@minimize, hyp, @gp, -20000, im, par1a{:}); 
    [L2 dL2] = feval(@gp,hyp_pN2, im, par1a{:});
    [ymu2 ys22 fmu2 fs22   ] = feval(@gp,hyp_pN2, im, par1b{:});    %Predict PGC (minus soma expression)
    
    %Now initialise hyperparameters over ESC    
    clear prior hyp
    hyp.cov    = [log(3);log(3)]; hyp.mean = mean(Y(:,1)); hyp.lik = log([4;4]);
    prior.mean = {[]}; prior.cov  = {pLS;[]}; prior.lik  = {[];[]};
    im = {@infPrior,@infVB,prior}; 
    par1a = {'meanConst',k1,'likT',X(find(X(:,2)==2),1),Y(find(X(:,2)==2))};
    par1b = {'meanConst',k1,'likT',X(find(X(:,2)==2),1),Y(find(X(:,2)==2)),Xstar(:,1)};    
    par1c = {'meanConst',k1,'likT',X(find(X(:,2)==2),1),Y(find(X(:,2)==2)),X(find(X(:,2)==2),1)}; %Predict at PGC    
    par1d = {'meanConst',k1,'likT',X(find(X(:,2)==2),1),Y(find(X(:,2)==2)),X(find(X(:,2)==3),1)}; %Predict at ESC locations        
    hyp_pN2b  = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         
    [L1b dL1b] = feval(@gp,hyp_pN1, im, par1a{:});
    [ymu1e  ys21e  fmu1e  fs21e] = feval(@gp,hyp_pN2b, im, par1b{:});
    [ymu1f  ys21f  fmu1f  fs21f] = feval(@gp,hyp_pN2b, im, par1c{:});
    [ymu1g  ys21g  fmu1g  fs21g] = feval(@gp,hyp_pN2b, im, par1d{:}); %This is predicting PGC at ESC locations
               
    clear prior hyp
    hyp.cov = [0.6;4;0.4;4;log(3);log(3)]; hyp.mean = []; hyp.lik = log([4;4]);
    prior.cov  = {cT3;pT;cT1;pT;pLS;[]}; prior.lik = {[];[]};
   
    im = {@infPrior,@infVB,prior}; 
    par1a = {'meanZero',k3,'likT',X(find(X(:,2)==3),1),Y(find(X(:,2)==3))-ymu1g};
    par1b = {'meanZero',k3,'likT',X(find(X(:,2)==3),1),Y(find(X(:,2)==3))-ymu1g,Xstar(:,1)};
    hyp_pN3 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});    
    [L3 dL3] = feval(@gp,hyp_pN3, im, par1a{:});     
    [ymu3 ys23 fmu3 fs23   ]= feval(@gp,hyp_pN3, im, par1b{:});  %Predict ESC (minus PGC expression)
    %%END INITIAL PARAMETER GUESSES
    
    
    %%Use ML to search over possible branching processes using initialised hyperparameters
    P_somapgc  = hyp_pN2.cov(1:2);
    P_pgcesc   = hyp_pN3.cov(1:4);                

    %Model 1: All bracnhes are different.    
    clear prior hyp    
    hyp.cov    = [P_somapgc; P_pgcesc; hyp_pN2.cov(3:4); hyp_pN3.cov(5:6); hyp_pN1.cov(1:2)];    
    hyp.mean   = mean(Y(:,1)); 
    hyp.lik    = max([hyp_pN1.lik,hyp_pN2.lik,hyp_pN3.lik]')';
    prior.mean = {[]};  
    prior.cov  = {cT2;pT;cT3;pT;cT1;pT;pLS;[];pLS;[];pLS;[]}; 
    prior.lik  = {[];[]};

    im         = {@infPrior,@infVB,prior};     
    par1a      = {'meanConst',k4,'likT',X,Y};
    par1b      = {'meanConst',k4,'likT',X,Y,Xstar};
    hyp_pN5    = feval(@minimize, hyp, @gp, -20000, im, par1a{:});
    [L5 dL5]   = feval(@gp,hyp_pN5, im, par1a{:});
    [ymu5 ys25 fmu5 fs25   ]= feval(@gp,hyp_pN5, im, par1b{:});
            
    
    %Model 2: Only soma is diff. 
    X1 = [x1,1*ones(size(x1,1),1); x2,2*ones(size(x2,1),1); x3,2*ones(size(x3,1),1)]; %Only soma different    
    hyp.cov  = [P_somapgc; [20;4;10;4]; hyp_pN2.cov(3:4); [5;5]; hyp_pN1.cov(1:2)];    
    hyp.mean = mean(Y(:,1)); 
    hyp.lik  = max([hyp_pN1.lik,hyp_pN2.lik,hyp_pN3.lik]')';
    prior.mean = {[]};  
    prior.cov  = {cT2;pT;pC;pC;pC;pC;pLS;[];pC;pC;pLS;[]}; 
    prior.lik  = {[];[]};    
 
    im       = {@infPrior,@infVB,prior};     
    par1a    = {'meanConst',k4,'likT',X1,Y};
    par1b    = {'meanConst',k4,'likT',X1,Y,Xstar};
    hyp_pN6  = feval(@minimize, hyp, @gp, -20000, im, par1a{:});  
    [L6 dL6] = feval(@gp,hyp_pN5, im, par1a{:});      
    [ymu6 ys26 fmu6 fs26   ]= feval(@gp,hyp_pN6, im, par1b{:});    

    %Model 3: Only ESCs are different  

    X3 = [x1,1*ones(size(x1,1),1); x2,1*ones(size(x2,1),1); x3,2*ones(size(x3,1),1)]; %Only soma different    
    hyp.cov  = [P_pgcesc; hyp_pN3.cov(5:6); hyp_pN1.cov(1:2)];    
    hyp.mean = mean(Y(:,1)); hyp.lik = max([hyp_pN1.lik,hyp_pN2.lik,hyp_pN3.lik]')';
    prior.mean = {[]};  
    prior.cov  = {cT3;pT;cT1;pT;pLS;[];pLS;[]}; 
    prior.lik  = {[];[]};
    
    im       = {@infPrior,@infVB,prior};     
    par1a    = {'meanConst',k5,'likT',X3,Y};
    par1b    = {'meanConst',k5,'likT',X3,Y,Xstar};
    hyp_pN8  = feval(@minimize, hyp, @gp, -20000, im, par1a{:});
    [L8 dL8] = feval(@gp,hyp_pN8, im, par1a{:});
    [ymu8 ys28 fmu8 fs28   ]= feval(@gp,hyp_pN8, im, par1b{:});
    
 
    
    %7) Nothing different     
    X6 = [x1,1*ones(size(x1,1),1); x2,1*ones(size(x2,1),1); x3,1*ones(size(x3,1),1)]; %Only soma different         

    hyp.cov    = [[20;4]; [20;4;10;4]; hyp_pN2.cov(3:4); hyp_pN3.cov(5:6); hyp_pN1.cov(1:2)];    
    hyp.cov    = [P_somapgc; P_pgcesc; hyp_pN2.cov(3:4); hyp_pN3.cov(5:6); hyp_pN1.cov(1:2)];    
    hyp.mean   = mean(Y(:,1)); 
    hyp.lik    = max([hyp_pN1.lik,hyp_pN2.lik,hyp_pN3.lik]')';
    prior.mean = {[]};  
    prior.cov  = {pC;pC;pC;pC;pC;pC;pC;pC;pC;pC;pLS;[]}; 
    prior.lik  = {[];[]};
 
    im         = {@infPrior,@infVB,prior};     
    par1a      = {'meanConst','covBranchingProcess_3Rec3','likT',X6,Y};
    par1b      = {'meanConst','covBranchingProcess_3Rec3','likT',X6,Y,Xstar};
    hyp_pN11   = feval(@minimize, hyp, @gp, -20000, im, par1a{:});
    [L11 dL11] = feval(@gp,hyp_pN11, im, par1a{:});
    [ymu11 ys211 fmu11 fs211   ]= feval(@gp,hyp_pN11, im, par1b{:});             

    keyboard
    %Store likelihoods
    L   = -[L1,L2,L1b,L3,L5,L6,L8,L11];
    BIC = -2*L + [4,6,4,8,14,8,10,4]*log(length(Y));
    ESCvED{i}.gene = genes{i};
    ESCvED{i}.L    = L;
    ESCvED{i}.BIC  = BIC;    
    
    ESCvED{i}.fmu1  = fmu1;
    ESCvED{i}.fs21  = fs21;    
    ESCvED{i}.ymu1  = ymu1;
    ESCvED{i}.ys21  = ys21; 
    
    ESCvED{i}.fmu2  = fmu2;
    ESCvED{i}.fs22  = fs22;    
    ESCvED{i}.ymu2  = ymu2;
    ESCvED{i}.ys22  = ys22;    

    ESCvED{i}.fmu3  = fmu3;
    ESCvED{i}.fs23  = fs23;    
    ESCvED{i}.ymu3  = ymu3;
    ESCvED{i}.ys23  = ys23;    
    
    %ESCvED{i}.fmu4  = fmu4;
    %ESCvED{i}.fs24  = fs24;    
    %ESCvED{i}.ymu4  = ymu4;
    %ESCvED{i}.ys24  = ys24;        

    ESCvED{i}.fmu5  = fmu5;
    ESCvED{i}.fs25  = fs25;    
    ESCvED{i}.ymu5  = ymu5;
    ESCvED{i}.ys25  = ys25;    
    
    ESCvED{i}.fmu6  = fmu6;
    ESCvED{i}.fs26  = fs26;    
    ESCvED{i}.ymu6  = ymu6;
    ESCvED{i}.ys26  = ys26;    

    %ESCvED{i}.fmu7 = fmu7;
    %ESCvED{i}.fs27 = fs27;    
    %ESCvED{i}.ymu7 = ymu7;
    %ESCvED{i}.ys27 = ys27;        

    ESCvED{i}.fmu8 = fmu8;
    ESCvED{i}.fs28 = fs28;    
    ESCvED{i}.ymu8 = ymu8;
    ESCvED{i}.ys28 = ys28;        
    
    %ESCvED{i}.fmu9 = fmu9;
    %ESCvED{i}.fs29 = fs29;    
    %ESCvED{i}.ymu9 = ymu9;
    %ESCvED{i}.ys29 = ys29;  
    
    %ESCvED{i}.fmu10 = fmu10;
    %ESCvED{i}.fs210 = fs210;    
    %ESCvED{i}.ymu10 = ymu10;
    %ESCvED{i}.ys210 = ys210;  
    
    ESCvED{i}.fmu11 = fmu11;
    ESCvED{i}.fs211 = fs211;    
    ESCvED{i}.ymu11 = ymu11;
    ESCvED{i}.ys211 = ys211;      
    
    ESCvED{i}.hyp1 = hyp_pN1;
    ESCvED{i}.hyp2 = hyp_pN2;    
    ESCvED{i}.hyp3 = hyp_pN3;    
    %ESCvED{i}.hyp4 = hyp_pN4;    
    ESCvED{i}.hyp5 = hyp_pN5;    
    ESCvED{i}.hyp6 = hyp_pN6;    
    %ESCvED{i}.hyp7 = hyp_pN7;    
    ESCvED{i}.hyp8 = hyp_pN8;    
    %ESCvED{i}.hyp9 = hyp_pN9;    
    %ESCvED{i}.hyp10 = hyp_pN10;    
    ESCvED{i}.hyp11 = hyp_pN11;               
    
    save(['ESCvED_7_' num2str(batchi) '_3Crtest.mat'],'ESCvED')    
    disp(['Step ' num2str(i)])
    catch
        ESCvED{i} = []; 
end

end

 save(['ESCvED_7_' num2str(batchi) '_3Crtest.mat'],'ESCvED') 

%keyboard
Output = 1;