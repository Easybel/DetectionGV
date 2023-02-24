%
%
% S2_GeneStat_wSim created by Isabel
%
%
%
%  What this script does: 
%  ---It takes the GeneOutput from the selection and no
%  selection condition and the GeneOut from the Nullmodel Script and plots
%  the multiple Hit statistic. 
% ---- It takes the function S2_CountMultiHits_fcn and finds the multi hits
% by its self
%
%%%%%%%%%%%%%%
%  Input: - GeneOutput for exp data
%         - GeneOut for sim data  
%         - bed file; masterlist
%       !!! You will want to exclude accessory genome genes (aka genes that are (partially) hit by acc. genome regions)
%         - Input: AccMM2Genes_Hits.mat - they are excluded if they have
%         the fraction exclude_thr
%
%        Additional input:
%             - you can exclude hits with a minimal fraction of minFrac
%%%%%%%%%%%%
%         Aditionally the script generates data to compare with from the very
%            simple binomial model
%%%%%%%%%%%%%%%%
%  Output
%         - all data is ordered into the ns and s structure
%         - figures are created where simulated and experimental, as well
%         as binomial model data, are compared
%

clear all; close all

 %% Load the paths 
A = load('/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/CNP2Genes/GeneOutput_Wns_cy20.mat');
ns.exp = A.GeneOutput;
A = load('/home/isabel/Dokumente/ExpEvol/Nullmodel/Wcy20/GeneOut_Wnsc20woRep2_3_10july2020.mat');
ns.sim = A.GeneOut;

A = load('/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/CNP2Genes/GeneOutput_Ws_cy20_allsamples.mat');
s.exp = A.GeneOutput;
A = load('/home/isabel/Dokumente/ExpEvol/Nullmodel/Wcy20/GeneOut_Wsc20woRep2_10_10july2020.mat');
s.sim = A.GeneOut;

% Other data
recipsize = 4215607;
masterlist = '/home/isabel/Dokumente/ExpEvol/dictionaries/BsubNC_000964wt/ml/BsubW23_2Bs166NCe_DP50_ml.txt';
recipbed = '/home/isabel/Dokumente/ExpEvol/dictionaries/BsubNC_000964wt/Bs166NCe_2020.bed.txt';
% deal with the acc genome
excAccmm = load('/home/isabel/Dokumente/ExpEvol/dictionaries/BsubNC_000964wt/acc/BsubW23_AccMM2Genes.mat');
excAccmm = excAccmm.AccMM2Genes;
exclude_thr = 1;

minFrac = 0.0001; % default is 0

%% Load data and handle it

% %Load recipient annonated bed file
fid = fopen(recipbed); bed = textscan(fid,'%s %f %f %s','delimiter',' ');
fclose(fid);
gene168.GN=bed{1}; gene168.S=bed{2}; gene168.E=bed{3}; gene168.BSU=bed{4}; 
gene168.L=gene168.E-gene168.S+1;
clear bed

%Load recipient/donor specific master list
fid = fopen(masterlist); imp = textscan(fid,'%f %s %s');
fclose(fid);
refchr.pos=imp{1}; refchr.atcg=imp{3}; 
clear imp
clear prompt str fid ans imp

idx_acc = [excAccmm([excAccmm(:).FracMean]>=exclude_thr).GeneNum];
mask_acc = ones(4422,1); mask_acc(idx_acc) = zeros(numel(idx_acc),1);
mask_acc = logical(mask_acc);
genes_woacc = sum(mask_acc);


%%  deal with the exp data and collect everything
% Do you want to know the MultiHitStat for all samples?
ns.dataset = 'Wns';
s.dataset = 'Ws';

ns.exp_samples = {'Wns0120','Wns0420','Wns0520','Wns0620','Wns0720','Wns0820'};
    
s.exp_samples = {'Ws0120','Ws0320','Ws0420','Ws0520','Ws0620','Ws0720','Ws0820', ...
           'Ws0920','Ws1120','Ws1220','Ws1320','Ws1420','Ws1520'};
       
ns.rep = numel(ns.exp_samples);
s.rep = numel(s.exp_samples);
       
% loop over the elements in Gs and only keep those that we are interested in
masks = zeros(numel([s.exp(:).RepNo]),1);
for i=1:numel([s.exp(:).RepNo])
    masks(i) = ~isempty(find(ismember({s.exp(i).Sample},s.exp_samples)));   
end
s.exp = s.exp(logical(masks)); clear masks
% loop over the elements in Gns and only keep those that we are interested in
maskn = zeros(numel([ns.exp(:).RepNo]),1);
for i=1:numel([ns.exp(:).RepNo])
    maskn(i) = ~isempty(find(ismember({ns.exp(i).Sample},ns.exp_samples)));   
end
ns.exp = ns.exp(logical(maskn)); clear maskn

% get the multihit data
ns.exp_multi = S2_CountMultiHits_fcn(excAccmm, ns.exp, gene168, refchr, exclude_thr);
s.exp_multi = S2_CountMultiHits_fcn(excAccmm, s.exp, gene168, refchr, exclude_thr);

% .. and get the number of genes hit per round
[a b] = unique([ns.exp(:).RepNo]); b = [b;numel([ns.exp(:).RepNo])+1];
ns.exp_Num = b(2:end) - b(1:end-1);
[c d] = unique([s.exp(:).RepNo]); d = [d;numel([s.exp(:).RepNo])+1];
s.exp_Num =  d(2:end) - d(1:end-1);

% and lastly exclude acc genes and those that are hit by a too small frac
ns.exp_woacc = ns.exp_multi(mask_acc);
s.exp_woacc = s.exp_multi(mask_acc);

%Here only the genes are taken into account that were hit at a fraction
%higher than minFrac
% all zeros are counted but a hit only counts if it at least minFrac of the gene
ns.exp_maskFrac = [ns.exp_woacc(:).FracMean]==0 | [ns.exp_woacc(:).FracMean]>minFrac;
ns.exp_SumHit = [ns.exp_woacc(ns.exp_maskFrac).SumHit];
ns.exp_HistV = histcounts(ns.exp_SumHit,-0.5:1:ns.rep+0.5);
ns.exp_SumSumHit = sum([ns.exp_SumHit]); % how many gene hit events are there?
% all zeros are counted but a hit only counts if it at least 1% of the gene
s.exp_maskFrac = [s.exp_woacc(:).FracMean]==0 | [s.exp_woacc(:).FracMean]>minFrac;
s.exp_SumHit = [s.exp_woacc(s.exp_maskFrac).SumHit];
s.exp_HistV = histcounts(s.exp_SumHit,-0.5:1:s.rep+0.5);
s.exp_SumSumHit = sum([s.exp_SumHit]); % how many gene hit events are there?

%%  deal with the sim data and collect everything
% Do you want to know the MultiHitStat for all samples?


for i=1:numel(ns.sim)
    ns.sim_multi(i).Out = S2_CountMultiHits_fcn(excAccmm, ns.sim(i).Out, gene168, refchr, exclude_thr);
    clear a b
    [a,b] = unique([ns.sim(i).Out.RepNo]); b = [b;numel([ns.sim(i).Out.RepNo])+1];
    ns.sim_Num(i,:) = b(2:end) - b(1:end-1);
    ns.sim_woacc(i).multi = ns.sim_multi(i).Out(mask_acc);
    
    ns.sim_maskFrac(i,:) = [ns.sim_woacc(i).multi.FracMean]==0 | [ns.sim_woacc(i).multi.FracMean]>minFrac;
    ns.sim_SumHit(i,:) = [ns.sim_woacc(i).multi(ns.sim_maskFrac(i,:)).SumHit];
    ns.sim_HistV(i,:) = histcounts(ns.sim_SumHit(i,:),-0.5:1:ns.rep+0.5);
    ns.sim_SumSumHit(i) = sum([ns.sim_SumHit(i,:)]); % how many gene hit events are there?

end
for i=1:numel(s.sim)
    s.sim_multi(i).Out = S2_CountMultiHits_fcn(excAccmm, s.sim(i).Out, gene168, refchr, exclude_thr);
    clear c d
    [c,d] = unique([s.sim(i).Out.RepNo]); d = [d;numel([s.sim(i).Out.RepNo])+1];
    s.sim_Num(i,:) =  d(2:end) - d(1:end-1);
    s.sim_woacc(i).multi = s.sim_multi(i).Out(mask_acc);
    
    s.sim_maskFrac(i,:) = [s.sim_woacc(i).multi.FracMean]==0 | [s.sim_woacc(i).multi.FracMean]>minFrac;
    s.sim_SumHit(i,:) = [s.sim_woacc(i).multi(s.sim_maskFrac(i,:)).SumHit];
    [s.sim_HistV(i,:) s.sim_Histx(i,:)] = histcounts(s.sim_SumHit(i,:),-0.5:1:s.rep+0.5);
    s.sim_SumSumHit(i) = sum([s.sim_SumHit(i,:)]); % how many gene hit events are there?

end

figure(1) % simulation data
ns.sim_Histerr1 =  std(ns.sim_HistV,1);
ns.sim_Histerr2 =  std(ns.sim_HistV,1);
h(1)=bar([0:ns.rep],mean(ns.sim_HistV,1),'FaceColor','b','EdgeColor','b','LineWidth',1,'FaceAlpha',0.2,'BarWidth',1);  hold on
er(1)=errorbar([0:ns.rep],mean(ns.sim_HistV,1),ns.sim_Histerr1,ns.sim_Histerr2,'LineStyle','none','Color','b')
h(2)=bar([0:ns.rep],ns.exp_HistV,'FaceColor','none','EdgeColor','r','LineWidth',1,'BarWidth',1);
legend([h(1) h(2)],{[ns.dataset ' sim'],[ns.dataset ' exp']}); 
ylim([0 1500]); %set(gca,'YScale','lin')
ylabel('Raw Hit Data'); xlabel('How often is the gene hit?')

figure(2) % simulation data
s.sim_Histerr1 =  std(s.sim_HistV,1);
s.sim_Histerr2 =  std(s.sim_HistV,1);
h(1)=bar([0:s.rep],mean(s.sim_HistV,1),'FaceColor','b','EdgeColor','b','LineWidth',1,'FaceAlpha',0.2,'BarWidth',1);  hold on
er(1)=errorbar([0:s.rep],mean(s.sim_HistV,1),s.sim_Histerr1,s.sim_Histerr2,'LineStyle','none','Color','b')
h(2)=bar([0:s.rep],s.exp_HistV,'FaceColor','none','EdgeColor','r','LineWidth',1,'BarWidth',1);
legend([h(1) h(2)],{[s.dataset ' sim'],[s.dataset ' exp']}); 
ylim([0 1400]); xlim([-1 12]); set(gca,'YScale','log')
ylabel('Raw Hit Data'); xlabel('How often is the gene hit?')

%% binomial model
allgenes = genes_woacc;
ns.binomial.ph = ns.exp_Num' / allgenes;  ns.binomial.pnh = (genes_woacc - ns.exp_Num') / allgenes;
m = 0;
for i=0:numel(ns.exp_Num) % this loop goes over all possible hit genes
m = m+1;
sett = 1:numel(ns.exp_Num);
idxHit = nchoosek([1:numel(ns.exp_Num)],i);
prop = 0;
for k = 1:size(idxHit,1)
idxnHit = setdiff(sett,idxHit(k,:));
if isempty(idxHit)
    prop = prop + prod([ns.binomial.pnh(idxnHit(k,:))]);
elseif isempty(idxnHit)
    prop = prop + prod([ns.binomial.ph(idxHit(k,:))]);
else 
    prop = prop + prod([ns.binomial.ph(idxHit(k,:)) ns.binomial.pnh(idxnHit)]);
end
end
ns.binomial.propHits(m) = i;
ns.binomial.prop(m) = prop;
end
ns.binomial.theory = round(ns.binomial.prop* allgenes);

s.binomial.ph = s.exp_Num' / allgenes;  s.binomial.pnh = (allgenes - s.exp_Num') / allgenes;
m = 0;
for i=0:numel(s.exp_Num) % this loop goes over all possible hit genes
m = m+1;
sett = 1:numel(s.exp_Num);
idxHit = nchoosek([1:numel(s.exp_Num)],i);
prop = 0;
for k = 1:size(idxHit,1)
idxnHit = setdiff(sett,idxHit(k,:));
if isempty(idxHit)
    prop = prop + prod([s.binomial.pnh(idxnHit(k,:))]);
elseif isempty(idxnHit)
    prop = prop + prod([s.binomial.ph(idxHit(k,:))]);
else 
    prop = prop + prod([s.binomial.ph(idxHit(k,:)) s.binomial.pnh(idxnHit)]);
end
end
s.binomial.propHits(m) = i;
s.binomial.prop(m) = prop;
end
s.binomial.theory = round(s.binomial.prop* allgenes);

figure(20)
b(1)=bar([0:ns.rep],mean(ns.sim_HistV,1),'FaceColor','b','EdgeColor','b','LineWidth',1,'FaceAlpha',0.2,'BarWidth',1);  hold on
ber(1)=errorbar([0:ns.rep],mean(ns.sim_HistV,1),ns.sim_Histerr1,ns.sim_Histerr2,'LineStyle','none','Color','b')
b(2)=bar([0:ns.rep],ns.binomial.theory,'FaceColor','none','EdgeColor','k','LineWidth',1,'FaceAlpha',0.7,'BarWidth',1);
b(3)=bar([0:ns.rep],ns.exp_HistV,'FaceColor','none','EdgeColor','r','LineWidth',1.5,'BarWidth',1);

legend([b(1) b(2) b(3)],{[ns.dataset ' sim'],[ns.dataset ' model prediction'],[ns.dataset ' exp data']}); 
ylim([0 1600]); xlim([-1 ns.rep+1]); set(gca,'YScale','lin')
ylabel('Hit Data'); xlabel('How often is the gene hit?')

figure(21)
b(1)=bar([0:s.rep],mean(s.sim_HistV,1),'FaceColor','b','EdgeColor','b','LineWidth',1,'FaceAlpha',0.2,'BarWidth',1);  hold on
ber(1)=errorbar([0:s.rep],mean(s.sim_HistV,1),s.sim_Histerr1,s.sim_Histerr2,'LineStyle','none','Color','b')
b(2)=bar([0:s.rep],s.binomial.theory,'FaceColor','none','EdgeColor','k','LineWidth',1,'FaceAlpha',0.7,'BarWidth',1);
b(3)=bar([0:s.rep],s.exp_HistV,'FaceColor','none','EdgeColor','r','LineWidth',1.5,'BarWidth',1);
legend([b(1) b(2) b(3)],{[s.dataset ' sim'],[s.dataset ' model prediction'],[s.dataset ' exp data']}); 
ylim([0 1600]); xlim([-1 s.rep]); set(gca,'YScale','log')
ylabel('Hit Data'); xlabel('How often is the gene hit?')
