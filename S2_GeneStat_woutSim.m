%
%
% S2_GeneStat_woutSim created by Isabel
%
%
%
%  What this script does: 
%  ---It takes the GeneOutput from the selection and no
%  selection condition and plots the multiple Hit statistic. 
% ---- It takes the function S2_CountMultiHits_fcn and finds the multi hits
% by its self
%
%%%%%%%%%%%%%%
%  Input: - GeneOutput for exp data
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
%         - figures are created where and experimental, as well
%         as binomial model data, are compared
%

clear all; close all

 %% Load the paths 
%A=  load('/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/CNP2Genes/GeneOutput_JWs_cy21.mat');

A = load('/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/CNP2Genes/GeneOutput_Wns_cy20.mat');
ns.exp = A.GeneOutput;
A = load('/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/CNP2Genes/GeneOutput_Ws_cy20_allsamples.mat');
s.exp = A.GeneOutput;

% Other data
recipsize = 4215607;
masterlist = '/home/isabel/Dokumente/ExpEvol/dictionaries/BsubNC_000964wt/ml/BsubW23_2Bs166NCe_DP50_ml.txt';
recipbed = '/home/isabel/Dokumente/ExpEvol/dictionaries/BsubNC_000964wt/Bs166NCe_2020.bed.txt';

 % deal with the acc genome
excAccmm = load('/home/isabel/Dokumente/ExpEvol/dictionaries/BsubNC_000964wt/acc/BsubW23_AccMM2Genes.mat');
excAccmm = excAccmm.AccMM2Genes;
exclude_thr = 1;

minFrac = 0.0001; % default is 0
method = 'linear'; % method for interpolation
clrmp = [0.1 0.1 0.1; 0 0.4 1; 0 0.79 0.61; 0.9 0 0; 0 0 0.6; 0 0 0.8; 0 0 1; 0 0.2 1; 0 0.4 1; 0 0.6 1]; 


%% Load data and handle it

%Load recipient annonated bed file
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


%% Do you want to know the MultiHitStat for all samples?
ns.dataset = 'Wns';
s.dataset = 'Ws';

ns.samples = {'Wns0120','Wns0420','Wns0520','Wns0620','Wns0720','Wns0820'};%,'Ws0720','Ws0820', ...
    %'Ws0920','Ws1020','Ws1120','Ws1220','Ws1320','Ws1420','Ws1520'};
    
s.samples = {'Ws0920','Ws1020','Ws1120','Ws1220','Ws1320','Ws1420'};%,'Ws1520'};

ns.rep = numel(ns.samples);
s.rep = numel(s.samples);

% loop over the elements in Gs and only keep those that we are interested in
masks = zeros(numel([s.exp(:).RepNo]),1);
for i=1:numel([s.exp(:).RepNo])
    masks(i) = ~isempty(find(ismember({s.exp(i).Sample},s.samples)));   
end
s.exp = s.exp(logical(masks)); clear masks
% loop over the elements in Gns and only keep those that we are interested in
maskn = zeros(numel([ns.exp(:).RepNo]),1);
for i=1:numel([ns.exp(:).RepNo])
    maskn(i) = ~isempty(find(ismember({ns.exp(i).Sample},ns.samples)));   
end
ns.exp = ns.exp(logical(maskn)); clear maskn

% get the multihit data
ns.multi = S2_CountMultiHits_fcn(excAccmm, ns.exp, gene168, refchr, exclude_thr);
s.multi = S2_CountMultiHits_fcn(excAccmm, s.exp, gene168, refchr, exclude_thr);

% .. and get the number of genes hit per round
[a b] = unique([ns.exp(:).RepNo]); b = [b;numel([ns.exp(:).RepNo])+1];
ns.Num = b(2:end) - b(1:end-1);
[c d] = unique([s.exp(:).RepNo]); d = [d;numel([s.exp(:).RepNo])+1];
s.Num =  d(2:end) - d(1:end-1);

% and lastly exclude acc genes and those that are hit by a too small frac
ns.woacc = ns.multi(mask_acc);
s.woacc = s.multi(mask_acc);

%Here only the genes are taken into account that were hit at a fraction
%higher than minFrac
% all zeros are counted but a hit only counts if it at least minFrac of the gene
ns.maskFrac = [ns.woacc(:).FracMean]==0 | [ns.woacc(:).FracMean]>minFrac;
s.maskFrac = [s.woacc(:).FracMean]==0 | [s.woacc(:).FracMean]>minFrac;

%% Now deal with the normalization and interpolation
% for ns
ns.SumHit = [ns.woacc(ns.maskFrac).SumHit];
% for all possible hits, find the bin values
ns.HistV = histcounts(ns.SumHit,-0.5:1:ns.rep+0.5);
ns.SumSumHit = sum([ns.SumHit]); % how many gene hit events are there?
% with which factor should the hits be normalized?
ns.normFactor = ns.SumSumHit;
ns.p = ns.SumHit / ns.normFactor;
% the interpolation starts at 1
ns.histcpx = [0:1:ns.rep]/ns.normFactor;
% ns.histcp = histc(ns.p,ns.histcpx);
% interpolate
ns.interp_px = linspace(0,max([ns.p]),100);
ns.interp_py = interp1(ns.histcpx(2:end),ns.HistV(2:end),ns.interp_px,method);

% and s
s.SumHit = [s.woacc(s.maskFrac).SumHit];
s.HistV = histcounts(s.SumHit,-0.5:1:s.rep+0.5);
s.SumSumHit = sum([s.SumHit]); % how many gene hit events are there?


% with which factor should the hits be normalized?
s.normFactor =  s.SumSumHit;
s.p = s.SumHit / s.normFactor;
% the interpolation starts at 1
s.histcpx = [0:1:s.rep]/s.normFactor;
%s.histcp = histc(s.p,s.histcpx);
% interpolate
s.interp_px = linspace(0,max([s.p]),100);
s.interp_py = interp1(s.histcpx(2:end),s.HistV(2:end),s.interp_px,method);


%% Plot multiple hit stat
figure(1)
h(2)=bar([1:6],s.HistV(2:end)./sum(s.Num),'FaceColor','b','EdgeColor','b','LineWidth',2,'FaceAlpha',0.2,'BarWidth',1);  hold on
h(1)=bar([1:6],ns.HistV(2:end)./sum(ns.Num),'FaceColor','none','EdgeColor','r','LineWidth',2,'BarWidth',1);
legend([h(1) h(2)],{ns.dataset,s.dataset}); 
ylim([0 1]); set(gca,'YScale','log')
ylabel('Raw Hit Data'); xlabel('How often is the gene hit?')

figure(11)
j(1)=plot(ns.HistV(1:end),ns.HistV(1:end),'ro','Displayname',[ns.dataset ' p_{rep}']);  hold on
j(2)=plot(s.HistV(1:end),s.HistV(1:end),'bo','Displayname',[s.dataset ' p_{rep}']);
ylim([1 max([ns.HistV s.HistV])+0.4*max([ns.HistV s.HistV])]); xlim([0 max([ns.histcpx s.histcpx])+0.1*max([ns.histcpx s.histcpx])]); 
    set(gca,'YScale','log');
j(3)=plot(ns.interp_px,ns.interp_py,':.','Color','r','Displayname','interpolated, lin');
j(4)=plot(s.interp_px,s.interp_py,':.','Color','b','Displayname','interpolated, lin');
ylabel('Hit Data'); xlabel('Normalized probability for a gene to be hit')
legend; 

%% binomial model
allgenes = genes_woacc;
ns.binomial.ph = ns.Num' / allgenes;  ns.binomial.pnh = (genes_woacc - ns.Num') / allgenes;
m = 0;
for i=0:numel(ns.Num) % this loop goes over all possible hit genes
m = m+1;
sett = 1:numel(ns.Num);
idxHit = nchoosek([1:numel(ns.Num)],i);
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

s.binomial.ph = s.Num' / allgenes;  s.binomial.pnh = (allgenes - s.Num') / allgenes;
m = 0;
for i=0:numel(s.Num) % this loop goes over all possible hit genes
m = m+1;
sett = 1:numel(s.Num);
idxHit = nchoosek([1:numel(s.Num)],i);
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

figure(20) %Monas idea
hold on
%title('Theory against data')
ns.fitp = polyfit(ns.binomial.theory,ns.HistV,1);
s.fitp=polyfit(s.binomial.theory,s.HistV,1);

ns.fitpwo0 = polyfit(ns.binomial.theory(ns.binomial.theory~=0),ns.HistV(ns.binomial.theory~=0),1);
s.fitpwo0 = polyfit(s.binomial.theory(s.binomial.theory~=0),s.HistV(s.binomial.theory~=0),1);
xvalues = linspace(0,2000,2000);

p(1)=plot(ns.binomial.theory,ns.HistV,'md','MarkerFaceColor','m','Displayname',ns.dataset');  hold on
pp(1) = plot(xvalues,xvalues*ns.fitp(1)+ns.fitp(2),'m--','Displayname',strcat('fit, m=',string(round(ns.fitp(1),3))))
p(2)=plot(s.binomial.theory,s.HistV,'bd','MarkerFaceColor','b','Displayname',s.dataset);
pp(2) = plot(xvalues,xvalues*s.fitp(1)+s.fitp(2),'b--','Displayname',strcat('fit, m=',string(round(s.fitp(1),3))))

xlabel('t_e: model prediction'); ylabel('hit data'); 
xlim([0 50]); ylim([0 75]); legend('Location','NorthWest'); 
plot([1 1e4],[1 1e4],'--','Color',[0.6 0.6 0.6],'Displayname','m = 1')


figure(21)
q(1)=bar([0:numel(ns.Num)],ns.HistV,'FaceColor','r','EdgeColor','r','LineWidth',2,'FaceAlpha',0.2,'BarWidth',1); hold on
q(2)=bar([0:numel(ns.Num)],ns.binomial.theory,'FaceColor','none','EdgeColor','k','LineWidth',2,'FaceAlpha',0.7,'BarWidth',1);
legend([q(1) q(2)],{[ns.dataset ' data'],[ns.dataset ' model prediction']}); 
ylim([0 3000]); xlim([-1 7]); set(gca,'YScale','log')
ylabel('Hit Data'); xlabel('How often is the gene hit?')

figure(22)
q(1)=bar([0:numel(s.Num)],s.HistV,'FaceColor','b','EdgeColor','b','LineWidth',2,'FaceAlpha',0.2,'BarWidth',1);   hold on
q(2)=bar([0:numel(s.Num)],s.binomial.theory,'FaceColor','none','EdgeColor','k','LineWidth',2,'FaceAlpha',0.2,'BarWidth',1);
legend([q(1) q(2)],{[s.dataset ' data'],[s.dataset ' model prediction']}); 
ylim([0 3000]); xlim([-1 13]); set(gca,'YScale','log')
ylabel('Hit Data'); xlabel('How often is the gene hit?')

figure(23) % Residual plot
q(1)=bar([0:numel(ns.Num)],ns.HistV-ns.binomial.theory,'FaceColor','none','EdgeColor','r','LineWidth',1,'FaceAlpha',0.2,'BarWidth',1);  
hold on
q(2)=bar([0:numel(s.Num)],s.HistV-s.binomial.theory,'FaceColor','b','EdgeColor','b','LineWidth',1,'FaceAlpha',0.2,'BarWidth',1);  hold on
legend([q(1) q(2)],{[ns.dataset ': data - model'],[s.dataset ': data - model']}); 
xlim([-1 13]); ylim([-220 250]); %ylim([-300 250]);
ylabel('Difference between measured and expected hits'); xlabel('How often is the gene hit?')

%% Look at identities
% Do identities correlate with p?

ns.ident = [ns.woacc(ns.maskFrac).mlGeneIdent];
s.ident = [s.woacc(s.maskFrac).mlGeneIdent];

% Correlation plot
figure(30)
scatter(ns.ident-0.93,ns.p,'xb','Displayname',ns.dataset)
hold on 
scatter(s.ident-0.93,s.p,'xr','Displayname',s.dataset)
legend('Location','NorthWest'); 
ylim([0 max([ns.p s.p])])
plot([0 0],[0 1],'k--','HandleVisibility','off'); xlabel('Difference to 93% Ident'); ylabel('per Gene, Hits/allGeneshit')
title('Ident vs. Hit probabilits per gene')
% Calculate correlation
ns.cova = cov(ns.ident,ns.p);
s.cova = cov(s.ident,s.p);
ns.corr = corrcoef(ns.ident,ns.p);
s.corr = corrcoef(s.ident,s.p);


%% Overview plots
% Plot Hits + Ident
% For this, again for all genes the ident and SumHit is needed.
ns.SumHitall = [ns.multi(:).SumHit];
ns.prepall = [ns.multi(:).SumHit]/ns.rep;
s.SumHitall = [s.multi(:).SumHit]; 
s.prepall = [s.multi(:).SumHit]/s.rep;
s.identall = [s.multi(:).mlGeneIdent];

totalgenes = 4422;
% Hot and cold plot - compare different exp. conditions
% be aware that the genes that appear here as acc. with ident=0 are the ones that 
% are excluded with the threshold
f12=figure('Units','centimeters','position',[20 20 16.97 24],'paperpositionmode','auto');   
rows = 6;
subplot(rows,1,1)
for i=1:rows
    figure(f12); subplot(rows,1,i)
    region = [floor(1+(i-1)*(totalgenes/rows)):ceil((totalgenes/rows)*i)];
    yyaxis left;
    plot(region,s.identall(region),'Color',[0.1 0.1 0.1 0.3],'HandleVisibility','off')
    ylim([0.8 1]); hold on
    yyaxis right
    if i==1
    bar(region,ns.prepall(region),'BarWidth',1,'FaceColor','r','displayname',ns.dataset)  
    hold on
    bar(region,s.prepall(region),'BarWidth',1,'FaceColor','b','displayname',s.dataset); legend('Location','NorthEastoutside')
    end
    bar(region,ns.prepall(region),'BarWidth',1,'FaceColor','r','HandleVisibility','off')
    hold on
    bar(region,s.prepall(region),'BarWidth',1,'FaceColor','b','HandleVisibility','off')
    xlim([1+(i-1)*(totalgenes/rows) (totalgenes/rows)*i ]); ylim([0 max([ns.prepall s.prepall])]);
    ax = gca; ax.YAxis(1).Color = [0.1 0.1 0.1]; ax.YAxis(2).Color = 'r';
end
subplot(rows,1,ceil(rows/2))
yyaxis left; ylabel('Identity')
yyaxis right; ylabel('p_{rep}')
subplot(rows,1,rows); xlabel('Gene number');

