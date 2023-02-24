%
%
%%%%%%%%%%%%%%%%   CNP2operons created by MonaIsa  %%%%%%%%%%%%%%%%%%%%%%%%
%
%  What this script does: It converts cluster information to operons,
%  counts how often operons are hit and does multipleHitStatistics for them
% 
%  Comment: Operons are named after their first BSUnumber (first on the genome)
%
%  Input: Cluster information; operon list; bed file (; masterlist)
%
%  Output
%  HitOutput: struct of every hit on a operon for every replicate
% -> HitOutput = struct('RepNo',[],'Cluster',[],'Type',[],'firstBSU',[],'firstGene',[],'opBSUs',[],'opGenes',[],'Frac',[]);  
%  Type: did the hit occur .. in, over a complete, at the front or
%         tail .. of a gene
%  OperonOutput: struct of every gene and how it was hit for every replicate
% -> OperonOutput = struct('RepNo',[],'Cluster',[],'Type',[],'firstBSU',[],'firstGene',[],'opBSUs',[],'opGenes',[],'Frac',[]);

%
%  RepSummary: summarizies information per replicate
% ->RepSummary = struct('RepNo',[],'Replicate',[],'OperonHit_fBSU',[],'OperonHit_fGene',[],'TotOperonHit',[]);
%
%  MultiHitStat = struct('Operon_fBSU',[],'Operon_fGene',[],'SumOperonHit',[],'IndexOperonOutput',[],'mlOperonIdent',[],'opBSUs',[]);

clear all; close all

% Define paths
masterlist = '/home/isabel/Dokumente/ExpEvol/scripts/CompareGenenames/W232Bsub.txt';
recipbed = '/home/isabel/Dokumente/ExpEvol/scripts/Genenames_Operons/BsubNC_000964wt_newRef_Mey.bed.txt';
operon = '/home/isabel/Dokumente/ExpEvol/scripts/Genenames_Operons/Operons_MeyIR_042020.mat';

% where are your replicates of interest in CNPSummary?  
dataset= 21:40;

% Path/names/suffix of the cluster cell:
CNPname = '/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/Cluster/CNPSummary_Wscy10_20_DP50.mat';
CNPSummary=load(CNPname);
CNPSummary=CNPSummary.CNPSummary;
CNPSummary = CNPSummary(dataset);

% the plaots and hotcolds are saved here
savepath = '/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/CNP2Genes/'

% %Load recipient/donor specific master list/ operon list
fid = fopen(masterlist); imp = textscan(fid,'%f %s');
fclose(fid);
refchr.pos=imp{1}; refchr.atcg=imp{3}; 
%Load recipient annonated bed file
fid = fopen(recipbed); bed = textscan(fid,'%s %f %f %f %f %s','delimiter',' ');
fclose(fid);
gene168I.GN=bed{1}; gene168I.S=bed{2}; gene168I.E=bed{3}; gene168I.L=bed{4};
gene168I.BSU=bed{6}; 
clear bed imp
clear prompt str fid ans
% %Load operon list

operon = load(operon);

prompt = 'Do you want to use operon_long [Y/N]';
str = input(prompt,'s');

if str == 'Y'
    operon = operon.Operons_filt_sort_long;
elseif str == 'N'
    operon = operon.Operons_filt_sort;
end

length_op = zeros(length(operon),1);
for i=1:length(operon)
    length_op(i,1) = length(operon(i).allBSUs);
end

prompt = 'How many genes should the operons have [value]';

x = input(prompt);

prompt = 'At least [A] or exact [E]';
str = input(prompt,'s');

if str == 'A'
    idx = find(length_op >= x);
elseif str == 'E'
    idx = find(length_op == x);
end
operon = operon(idx);


%User input: Use which cluster information?
prompt = 'Do you want to use C or Adist to find Cluster? C/Adist [C]: ';
str = input(prompt,'s');

%User input: Use which cluster information?

prompt = 'Do you want to use C or Adist to find Cluster? C/Adist [C]: ';
which = input(prompt,'s');
m=0;
if which=='C'
    C={CNPSummary(:).C}';
    for i=1:size(C,1)
        maskSample(i) = ~isempty(C{i,1});
        if ~isempty(C{i,1})
            m = m+1;
            Cluster{m,1}(1,:) = C{i,1}(1,:);
            Cluster{m,1}(2,:) = C{i,1}(2,:);
        end
    end
    
elseif isempty(which) || true(which == 'Adist')
    Adist={CNPSummary(:).Adist}';
    for i=1:size(Adist,1)
        maskSample(i) = ~isempty(Adist{i,1});
        if ~isempty(Adist{i,1})
            m = m+1;
            Cluster{m,1}(1,:) = Adist{i,1}(:,2)';
            Cluster{m,1}(2,:) = Adist{i,1}(:,3)';
        end
    end
end

%%%%%Load variables%%%%%%%%%%%%
recipsize = 4215607;
% Here write the names of the replicates that the clusters refer to 
rep = {CNPSummary(maskSample).Samples};
ORI = {CNPSummary(maskSample).ORI_Crossing};

%% HitOutput is generated 
%Here for every replicate i all cluster c are checked for hits with operons and all different hits are written as entry into HitOutput.
%   -> This means: operons can appear multiple times per replicate -> this is
%      fixed in operonOutput
HitOutput = struct('RepNo',[],'Cluster',[],'Type',[],'firstBSU',[],'firstGene',[],'opBSUs',[],'opGenes',[],'Frac',[]);
m = 0;
for i=1:size(Cluster,1)    
    % if we have a ORI crossing, split up the first entry in 2:
    % start(1):recipsize & 1:end(1)
    if ORI{i} == 1
        first = [1 Cluster{i}(2,1)]; last = [Cluster{i}(1,1) recipsize]; 
        Cluster{i}(:,1) = first; Cluster{i}(:,end+1) = last;
    end
%   
    for c=1:size(Cluster{i},2)        
        ptail = find( ((Cluster{i}(1,c) > [operon(:).S]) & (Cluster{i}(1,c) < [operon(:).E]) & (Cluster{i}(2,c) >= [operon(:).E])) == 1);      
        for pt=1:numel(ptail)
            m = m+1;
            HitOutput(m).Type = 'tail';
            HitOutput(m).RepNo = i;   HitOutput(m).Cluster = c;
            HitOutput(m).firstBSU = {operon(ptail(pt)).firstBSU}; HitOutput(m).firstGene = operon(ptail(pt)).allGenes(1); 
            HitOutput(m).opBSUs = operon(ptail(pt)).allBSUs; HitOutput(m).opGenes = operon(ptail(pt)).allGenes; 
            HitOutput(m).Frac = ( operon(ptail(pt)).E - Cluster{i}(1,c) + 1)/operon(ptail(pt)).L;          
        end        
        complete = find( ((Cluster{i}(1,c) <= [operon(:).S]) & (Cluster{i}(2,c) >= [operon(:).E])) == 1 );
        for co=1:numel(complete)
            m = m+1;
            HitOutput(m).Type = 'complete';
            HitOutput(m).RepNo = i;    HitOutput(m).Cluster = c;
            HitOutput(m).firstBSU = {operon(complete(co)).firstBSU}; HitOutput(m).firstGene = operon(complete(co)).allGenes(1); 
            HitOutput(m).opBSUs = operon(complete(co)).allBSUs; HitOutput(m).opGenes = operon(complete(co)).allGenes; 
            HitOutput(m).Frac = 1;
        end       
        pin = find(  ((Cluster{i}(1,c) > [operon(:).S]) & (Cluster{i}(2,c) < [operon(:).E]))==1  );
        for pi=1:numel(pin)
            m = m+1;
            HitOutput(m).Type = 'in';
            HitOutput(m).RepNo = i;    HitOutput(m).Cluster = c;
            HitOutput(m).firstBSU = {operon(pin(pi)).firstBSU}; HitOutput(m).firstGene = operon(pin(pi)).allGenes(1); 
            HitOutput(m).opBSUs = operon(pin(pi)).allBSUs; HitOutput(m).opGenes = operon(pin(pi)).allGenes; 
            HitOutput(m).Frac = (Cluster{i}(2,c) - Cluster{i}(1,c) + 1)/operon(pin(pi)).L;
        end       
        pfront = find( ((Cluster{i}(1,c) <= [operon(:).S]) & (Cluster{i}(2,c) >  [operon(:).S]) & (Cluster{i}(2,c) < [operon(:).E]))==1);      
        for pf=1:numel(pfront)
            m= m+1;
            HitOutput(m).Type = 'front';
            HitOutput(m).RepNo = i;    HitOutput(m).Cluster = c;
            HitOutput(m).firstBSU = {operon(pfront(pf)).firstBSU}; HitOutput(m).firstGene = operon(pfront(pf)).allGenes(1); 
             HitOutput(m).opBSUs = operon(pfront(pf)).allBSUs; HitOutput(m).opGenes = operon(pfront(pf)).allGenes;
            HitOutput(m).Frac = (Cluster{i}(2,c) - operon(pfront(pf)).S + 1)/operon(pfront(pf)).L;      
        end       
        clear complete pin pfront ptail pi pf pt co 
    end
end

%% GeneOutput and RepSummary are generated
% Up to now in OerponOutput, Operons can be hit multiple times in a certain
% replicate with
% what has to be considered: only BSU numbers are unique! genenames are
% not, so this here is based on BSU numbers
OperonOutput = struct('RepNo',[],'Cluster',[],'Type',[],'firstBSU',[],'firstGene',[],'opBSUs',[],'opGenes',[],'Frac',[]);
RepSummary = struct('RepNo',[],'Replicate',[],'OperonHit_fBSU',[],'OperonHit_fGene',[],'TotOperonHit',[]);

l = 0;
for i=1:numel(rep)
idx = find([HitOutput.RepNo] == i);
b = {HitOutput(idx).firstBSU}'; BSUcollectRep = vertcat(b{:});
for j=1:numel(idx)
match = find(strcmp(BSUcollectRep{j},BSUcollectRep));
matchD = match(match~=j);
if length(match)==1
l = l + 1;
OperonOutput(l) = HitOutput(idx(j));
elseif matchD>j
l = l + 1;
OperonOutput(l) = HitOutput(idx(j));
%add the other gene infos to this entry
OperonOutput(l).Frac = sum([HitOutput(idx(match)).Frac]);
OperonOutput(l).Cluster = [HitOutput(idx(match)).Cluster];
OperonOutput(l).Type = {HitOutput(idx(match)).Type};
end
end

%Now for each Replicate I collect the summary
idxUniq = find([OperonOutput.RepNo] == i);
bb = {OperonOutput(idxUniq).firstBSU}'; RepSummary(i).OperonHit_fBSU = vertcat(bb{:});
bbb = {OperonOutput(idxUniq).firstGene}'; RepSummary(i).OperonHit_fGene = vertcat(bbb{:});
RepSummary(i).RepNo = i; RepSummary(i).Replicate = rep{i};
RepSummary(i).TotOperonHit = length(idxUniq);
end
clear b bb bbb match matchD idxUniq
%% Multiple HitStatistics: How often is each Gene hit over all replicates??

MultiHitStat = struct('Operon_fBSU',[],'Operon_fGene',[],'SumOperonHit',[],'IndexOperonOutput',[],'mlOperonIdent',[],'opBSUs',[]);
c = {OperonOutput(:).firstBSU}';          BSUcollect = vertcat(c{:});
for i=1:numel([operon(:).L])
    MultiHitStat(i).Operon_fBSU = operon(i).firstBSU;
    MultiHitStat(i).Operon_fGene = operon(i).allGenes(1);
    MultiHitStat(i).opBSUs = operon(i).allGenes;
    idx = find(strcmp(OperonOutput(i).firstBSU,BSUcollect));
    if ~isempty(idx)
        MultiHitStat(i).SumOperonHit = numel(idx);
        MultiHitStat(i).IndexOperonOutput = idx;
        clear idx
    elseif isempty(idx)
        MultiHitStat(i).SumOperonHit = 0;
    end
    
    %add the gene Identity based on the masterlist to MultiHitStat
    MultiHitStat(i).mlOperonIdent =  1-(sum(ismember(refchr.pos,[operon(i).S:operon(i).E]))/operon(i).L);
end
clear c 

%% Hot and cold plot
rownum = 6; operonnum = size([operon(:).L],2); maxMult = max([MultiHitStat(:).SumOperonHit]);
figure(1)
for i=1:rownum
    subplot(6,1,i)
    plot([1:operonnum],[MultiHitStat(:).mlOperonIdent]','r')
    hold on
    yyaxis right
    bar([1:operonnum],[MultiHitStat(:).SumOperonHit])
    hold on
    xlim([1+(i-1)*(operonnum/rownum) (operonnum/rownum)*i ])
    ylim([0 maxMult+2])
end
subplot(6,1,1)
title('Hot and Cold plot - Operons')
%% Plot operon positions
rownum = 10; operonnum = size([operon(:).L],2);
figure(2)
for i=1:rownum
subplot(10,1,i)
    for j=1:operonnum
        if all(operon(j).S > [operon(1:j-1).E])
        m = 1;
        elseif sum(operon(j).S < [operon(1:j-1).E])>0
            m = m + 1;
        end
        line([operon(j).S operon(j).E],[m m],'LineWidth',3)
        
        hold on    
        xlim([1+(i-1)*(recipsize/rownum) (recipsize/rownum)*i ])
    ylim([0 6])
    end
    
end
subplot(10,1,1)
title('Operons across the genome')
%% Multi Hit Statistics Plots
figure(3)
histogram([MultiHitStat(:).SumOperonHit])
set(gca, 'YScale', 'lin')
title('Multiple Hit Statistic')
ylim([0 500])

% figure(4)
% histogram([MultiHitStat(:).SumOperonHit]/(operonnum*numel(rep)))
% set(gca, 'YScale', 'log')


%% Save structures
% save('RepSummary_Jeff2BsNCe_Adist_oldmlacc.mat','RepSummary')
% save('Cluster_Jeff2BsNCe_Adist_oldmlacc.mat','Cluster')
% save('HitOutput_Jeff2BsNCe_Adist_oldmlacc.mat','HitOutput')
% save('GeneOutput_Jeff2BsNCe_Adist_oldmlacc.mat','GeneOutput')
% 
