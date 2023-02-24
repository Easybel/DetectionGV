function [MultiHitStat] = S2_CountMultiHits_fcn(excAccmm, GeneOutput, gene168, refchr, exclude_thr)
%% Multiple HitStatistics: How often is each Gene hit over all replicates??
% Prepare gene list so that it excludes the genes that are excAccmm and 
%that are above the threshold
acm = {{excAccmm(:).BSU}'}; accmm.BSU = vertcat(acm{:});

mask = [excAccmm(:).FracMean]>=exclude_thr;
ex = {{excAccmm(mask).BSU}'}; exc.BSU = vertcat(ex{:});

mlSNPsinGene = 0;
MultiHitStat = struct('BSU',[],'Genename',[],'SumHit',[],'HitIndex',[],'FracMean',[],'mlGeneIdent',[],'Samples',[]);
c = {GeneOutput(:).BSU}';          BSUcollect = vertcat(c{:});
Samplecollect = {GeneOutput(:).Sample}';   
n = 0;
% here, i goes through all the entries in the bed file
for i=1:numel(gene168.BSU)
    
    MultiHitStat(i).BSU = gene168.BSU{i};
    MultiHitStat(i).Genename = gene168.GN{i};
    MultiHitStat(i).GeneNum = i;
    %add the gene Identity of the i-th gene based on the masterlist to
    %MultiHitStat- if the i-th gene also appears on the excAnnMM list, then
    %calculate it a bit differently
    idx_acc = find(strcmp(gene168.BSU{i},accmm.BSU));
    if ~isempty(idx_acc)
        % here a very small number is added just to not to devide by 0
        l = gene168.L(i)*(1-(excAccmm(idx_acc).FracMean));
        snps = sum(ismember(refchr.pos,gene168.S(i):gene168.E(i)));
        MultiHitStat(i).mlGeneIdent =  1-(snps/l);
        if snps>l; snps = l; MultiHitStat(i).mlGeneIdent =  1-(snps/l); end
        if l == 0; l=0.0000001; MultiHitStat(i).mlGeneIdent =  0; end        
    else
        MultiHitStat(i).mlGeneIdent =  1-(sum(ismember(refchr.pos,gene168.S(i):gene168.E(i)))/gene168.L(i));
    end
    
% if the gene also appears in exc.BSU, mark it as excluded in the
% accxluded field
    exc_idx = find(strcmp(gene168.BSU{i},exc.BSU));
% at which index idx did the i-th gene fit the geneoutput list?
    idx = find(strcmp(gene168.BSU{i},BSUcollect));
    if ~isempty(exc_idx)
        MultiHitStat(i).accxluded = 1;
        MultiHitStat(i).SumHit = 0;
        MultiHitStat(i).FracMean = 0;
        
% if the gene was hit at least once, calculate ..
    elseif isempty(exc_idx) && ~isempty(idx)
        MultiHitStat(i).accxluded = 0;
        MultiHitStat(i).SumHit = numel(idx);
        MultiHitStat(i).HitIndex = idx;
        MultiHitStat(i).FracList = [GeneOutput(idx).Frac];
        MultiHitStat(i).FracMean = mean([GeneOutput(idx).Frac]);
        MultiHitStat(i).Samples = {Samplecollect{idx}};
        
    elseif isempty(exc_idx) && isempty(idx)
        MultiHitStat(i).accxluded = 0;
        MultiHitStat(i).SumHit = 0;
        MultiHitStat(i).FracMean = 0;
    end
    
    clear idx
end
end