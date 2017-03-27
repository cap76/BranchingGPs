function [Output] = Stella_allgenes(Batchi)

%Load the pseudotimes
load(['../results/primordial_germ_cells/Pseudotime/Stella_Pseudotime_' num2str(1) '__2E_final.mat'])                 

%Now load in DE genes
%DE = importdata('StellaDE.csv',',',1);
DE = importdata('./DEall.csv',',',0)

%Now load in the expression data
%D1 = importdata('Mouse_Yun_1a.txt');
%D2 = importdata('Mouse_Yun_1b.txt');
%D1.data(2:end,:) = log2(D1.data(2:end,:)+1);
%D2.data(2:end,:) = log2(D2.data(2:end,:)+1);  
%D1.data = [D1.data,D2.data]';
D1 = importdata('./Figure_1.xlsx')
D1.data(3:end,:) = log2(D1.data(3:end,:)+1);

%D1 = importdata('GSE80810_RPRT_WT.csv');
%D2 = importdata('GSE80810_RPRT_KO.csv');
%D1.data(2:end,:) = log2(D1.data(2:end,:)+1);
%D2.data(2:end,:) = log2(D2.data(2:end,:)+1);
% Data is in the following format:
% (No. genes + 3) x (No. cells) 
%
% Row 1 (Type): ESC/unlabelled (-1) Preimplantation (0) PGC (1) or Soma (2)
% Row 2 (Sex):  Unkown (0) Female (1) Male (2)
% Row 3 (Capture time): ESC (-1), Preimplantation 0-6, specfied cell PGC or soma (>6)
% Row 4 onwards: Gene expression levels
inds0 = find(D1.data(1,:)==0 & D1.data(2,:)==1);
inds1 = find(D1.data(1,:)==1 & D1.data(2,:)==1);
inds2 = find(D1.data(1,:)==3 & D1.data(2,:)==1);
inds3 = find(D1.data(1,:)==0 & D1.data(2,:)==2);
inds4 = find(D1.data(1,:)==1 & D1.data(2,:)==2);    
inds5 = find(D1.data(1,:)==3 & D1.data(2,:)==2);       
t1   = [D1.data(1,[inds0,inds1,inds2]),D1.data(1,[inds3,inds4,inds5])];            %Vector of developmental stage (capture time)
D1.data = D1.data(:,[inds0,inds1,inds2,inds3,inds4,inds5]);
uID  = 1:1:length(t1);          %Unique ID (for rejumbling everything later on)


%Get X,Y in correct format
Ya     = D1.data(3:end,:)';
X     = Output.Param.Store{end}.update.x;

%Prediction points
Xstar = [linspace(min(X(:,1)),max(X(:,1)),1000)',ones(1000,1);linspace(min(X(:,1)),max(X(:,1)),1000)',ones(1000,1)];

%Covariance functions
k1 = {'covMaterniso',3};
k2 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}};

startin = (Batchi-1)*500+1;
endin = min(Batchi*1000,28155)

for arc = startin:endin%1:size(Ya,2) %Loop over all DE genes
   
    %Find the DE score of 
    gene = D1.textdata{arc+3,1};
    ind  = find(strcmp(gene,DE.textdata(2:end,1))==1);
    DEscore = min(DE.data(ind,end-2:end));
    
    Y = Ya(:,arc);
    
    if DEscore<0.5
    
        try
    NZI(arc,1) = 1;
        
    %Some prior distributions 
    cT  = {@priorGauss,0.1,1}; %CP
    pLS = {@priorGauss,0.5,1}; %Prior over length scale.
    pN  = {@priorGauss,0.5,1}; %Noise
    pT  = {@priorGauss,3,1};   %Branch rate
    
    %First fit base process (WT) to guess some params
    hyp.cov = [log(3);log(3)]; hyp.mean = mean(Y); hyp.lik = log(2);
    prior.mean = {[]};  prior.cov  = {pLS;[]}; prior.lik = {pN};
    im    = {@infPrior,@infExact,prior}; 
    par1a = {'meanConst',k1,'likGauss',X(find(X(:,2)==1),1),Y(find(X(:,2)==1))};
    par1b = {'meanConst',k1,'likGauss',X(find(X(:,2)==1),1),Y(find(X(:,2)==1)),Xstar(:,1)};
    par1c = {'meanConst',k1,'likGauss',X(find(X(:,2)==1),1),Y(find(X(:,2)==1)),X(find(X(:,2)==2),1)}; %Predict here    
    hyp_pN1 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         
    [L1 dL1] = feval(@gp,hyp_pN1, im, par1a{:});         
    [ymu1 ys21 fmu1 fs21   ]= feval(@gp,hyp_pN1, im, par1b{:});
    [ymu1c ys21c fmu1c fs21c   ]= feval(@gp,hyp_pN1, im, par1c{:});

    %Now guess Stella params.
    clear prior hyp
    hyp.cov = [0.4;4;log(3);log(3)]; hyp.lik = log(2);
    prior.cov  = {cT;pT;pLS;[]}; prior.lik = {pN};
    im = {@infPrior,@infExact,prior}; 
    par1a = {'meanZero',k2,'likGauss',X(find(X(:,2)==2),1),Y(find(X(:,2)==2))-ymu1c};
    par1b = {'meanZero',k2,'likGauss',X(find(X(:,2)==2),1),Y(find(X(:,2)==2))-ymu1c,Xstar(:,1)};
    hyp_pN2 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L2 dL2] = feval(@gp,hyp_pN2, im, par1a{:});         % optimise    
    [ymu2 ys22 fmu2 fs22   ]= feval(@gp,hyp_pN2, im, par1b{:});    
    
    %Now do full optimsation of the branch process
    clear prior hyp
    hyp.cov = [hyp_pN2.cov(1:4);hyp_pN1.cov(1:2)]; hyp.mean = hyp_pN1.mean; hyp.lik = log(2);
    prior.cov  = {cT;pT;pLS;[];pLS;[]}; prior.lik = {pN}; prior.mean = {[]};
    im = {@infPrior,@infExact,prior}; 
    par1a = {'meanConst','covBranchingProcess_2E','likGauss',X,Y};
    par1b = {'meanConst','covBranchingProcess_2E','likGauss',X,Y,Xstar};
    hyp_pN3 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L3 dL3] = feval(@gp,hyp_pN3, im, par1a{:});         % optimise    
    [ymu3 ys23 fmu3 fs23   ]= feval(@gp,hyp_pN3, im, par1b{:});    
    
    
    %Now do full optimsation of joint process
    clear prior hyp
    hyp.cov = [hyp_pN1.cov(1:2)]; hyp.mean = hyp_pN1.mean; hyp.lik = hyp_pN1.lik;
    prior.cov = {pLS;[]}; prior.lik = {pN}; prior.mean = {[]};
    im = {@infPrior,@infExact,prior}; 
    par1a = {'meanConst',k1,'likGauss',X(:,1),Y};
    par1b = {'meanConst',k1,'likGauss',X(:,1),Y,Xstar(:,1)};
    hyp_pN4 = feval(@minimize, hyp, @gp, -20000, im, par1a{:});         % optimise
    [L4 dL4] = feval(@gp,hyp_pN4, im, par1a{:});         % optimise    
    [ymu4 ys24 fmu4 fs24   ]= feval(@gp,hyp_pN4, im, par1b{:});        
    
    L = -[L1,L2,L3,L4];
    Stella{arc}.gene = gene;
    Stella{arc}.L    = L;
    Stella{arc}.BIC  = -2*L + [4,6,8,4]*log(66);

    Stella{arc}.H1 = hyp_pN1;
    Stella{arc}.H2 = hyp_pN2;
    Stella{arc}.H3 = hyp_pN3;    
    Stella{arc}.H4 = hyp_pN4;    

    
    Stella{arc}.fmu1 = fmu1;
    Stella{arc}.fs21 = fs21;    
    Stella{arc}.ymu1 = ymu1;
    Stella{arc}.ys21 = ys21;    
    Stella{arc}.fmu2 = fmu2;
    Stella{arc}.fs22 = fs22;    
    Stella{arc}.ymu2 = ymu2;
    Stella{arc}.ys22 = ys22;
    Stella{arc}.fmu3 = fmu3;
    Stella{arc}.fs23 = fs23;    
    Stella{arc}.ymu3 = ymu3;
    Stella{arc}.ys23 = ys23;  
    Stella{arc}.fmu4 = fmu4;
    Stella{arc}.fs24 = fs24;    
    Stella{arc}.ymu4 = ymu4;
    Stella{arc}.ys24 = ys24;      
    
     disp(['Step ' num2str(arc) ' of ' num2str(size(Ya,2))])
    
    save(['Stella_AllBranching_Matern_re'  num2str(Batchi) '.mat'],'Stella')   
        catch %Numerically unsound
        NZI(arc,1) = 0;
        Stella{arc} = [];
        end
    else %Not DE
    NZI(arc,1) = 0;
    Stella{arc} = [];
    end
end

%inds = find(NZI==1);
%Stella = Stella(inds);

save(['Stella_AllBranching_Matern_re' num2str(Batchi) '_c.mat'],'Stella')