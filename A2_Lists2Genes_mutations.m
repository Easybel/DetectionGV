%
%
%%%%%%%%%%%%%%%%   CNP2genes created by MonaIsa  %%%%%%%%%%%%%%%%%%%%%%%%
%
%  What this script does: It converts cluster information to genenames,
%  counts how often genes are hit and does multipleHitStatistics for them
%
%%%%%%%%%%%%%%
%  Input: - Cluster information: either CNPSummary, deldup or a txt file
%         - bed file; masterlist
%         - dataset: where in SNPSummary or deldup are the data of
%         interest?
%         - savepath: where are the data saved?
%       !!! You will want to exclude accessory genome genes (aka genes that are (partially) hit by acc. genome regions)
%         - Input: AccMM2Genes_Hits.mat
%
%        Additional input:
%             - you can load a txt file with a lsit of BSU names of genes that
%               interest you and do the whole analysis just for these genes, otherwise
%               it is done with all genes
%
%%%%%%%%%%%%%%%%
%  Output
%         HitOutput: struct of every hit on a gene for every replicate
%             -> struct('RepNo',[],'Cluster',[],'Type',[],'BSU',[],'Genename',[],'Frac',[],'IdentMl',[]);
%         Type: did the hit occur .. in, over a complete, at the front or
%         tail .. of a gene
%
%          GeneOutput: struct of every gene and how it was hit for every replicate
%             -> struct('RepNo',[],'Cluster',[],'Type',[],'BSU',[],'Genename',[],'Frac',[],'IdentMl',[]);
%
%          RepSummary: summarizies information per replicate
%              -> struct('RepNo',[],'BSUHit',[],'GeneHit',[],'Replicate',[],'UniqBSU',[]);
%
%          MultiHitStat = struct('BSU',[],'Genename',[],'SumHit',[],'HitIndex',[]);
%
%
% things to keep in mind:
%    - the genes that belong to accessory genome have to be excluded from
%    further analysis. For this reason

clear all; close all

%% Here, you have to make some inputs:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Data Input   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% different data inputs, choose one:
IndataType = "txtListStartEnd";       % or "Adist" or "denovo" --> they all come from CNPSummary
                            % or "dup" "del"         --> from Output_DelDup
                            % or "txtListStart" for a list of mutations e.g.
                            % or "txtListStartEnd"

% Here give the data file
% either: CNPSummary.mat, Output_deldup.mat or txt document

IndataPath  = "/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/MS11_ncbi_IR/MULTIMAPPER_Ms11toMs11/";
%"/home/isabel/Dokumente/P5_ExpEvol_Ngo/DNA/analysis/Mutations_2_Genes/";
%"/home/isabel/Dokumente/P2_FitnessExp/0_Libraries/1_LibDataSeq_Bval/2_Analysis/";
%"/home/isabel/Dokumente/P5_ExpEvol_Ngo/DNA/data/";
%"/home/isabel/Dokumente/P2_FitnessExp/0_Libraries/1_LibDataSeq_Bval/2_Analysis/";
Indata      = "Blast_MS11Genes_2_MS11_editted.tab";
%"mlNgoDG4_2_MS11.txt";
%"LibSCBval_filt_20210121_CNPSummary.mat";
%"mlNgoDG4_2_MS11.txt";
%"LibSCBval_filt_20210121_CNPSummary.mat";
%dataInput = '/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/Del_Dup/Output_deldup_Wscy20_woRep10_15julyIR.mat';

% If you hand a .mat, where are your replicates of interest?
dataset= [1:5]; %[36:40];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Bed File     %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now we need the information, where which gene is in your organism
recipbed   = "/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/MS11_ncbi_IR/Manual_Lists/MS11.bed.mat";
%"/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/0_dictionaries/Bs166_unchanged/NC000964_NCBI_unchanged.bed.mat";
%"/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/MS11_ncbi_IR/Manual_Lists/MS11.bed.mat";
masterlist = '/home/isabel/Dokumente/ExpEvol/dictionaries/BsubNC_000964wt/ml/BsubW23_2Bs166NCe_DP50_ml.txt';
%Path to the MultiHitStat output wth the acc and mm genes
% excAccmm = load('/home/isabel/Dokumente/ExpEvol/dictionaries/BsubNC_000964wt/acc/BsubW23_AccMM2Genes.mat');
% excAccmm = excAccmm.AccMM2Genes;
% exclude_thr = 1.1;

% Do you want to only search for a subset of the genes given in the
% bed?
SearchInSpecGenes = "OFF";
specGenes = "/home/isabel/Dokumente/P5_ExpEvol_Ngo/DNA/analysis/RandomGenes.txt";

%%%%%Load variables%%%%%%%%%%%%
recipsize = 4215607;
yes=1; no=0;
% the plots and hotcolds are saved here
savepath = '/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/denovos2Genes/'
%% Here the data is loaded

% %Load recipient annonated bed file
bed     = load(recipbed);
bedFile_tmp = bed.uniqList;


if SearchInSpecGenes == "ON"
    % load special genes
    fid = fopen(specGenes);
    imp = textscan(fid,'%s','delimiter','\t');
    fclose(fid);
    specGen.N=imp{1};
    clear imp fid
    
    % and look for the indices of the special genes
    [Lia,Locb] = ismember([bedFile_tmp(:).locus_tag],{specGen.N{:}});
    speIdx = find(Lia);
    bedFile = bedFile_tmp(speIdx)
else
    bedFile = bedFile_tmp;
end
tmp = struct("L", num2cell([bedFile.ende] - [bedFile.start] +1));
[bedFile.L] = deal(tmp.L)
% Read in data .. depending on dataType

if ismember(IndataType,["Cdist","Adist","denovo"])
    CNPSummary=load(IndataPath + Indata);
    CNPSummary=CNPSummary.CNPSummary(dataset);
    C     = {CNPSummary(:).C}';
    Adist = {CNPSummary(:).Adist}';
    denovo= {CNPSummary(:).denovo}';
    if IndataType == "Cdist"
        for i=1:size(C,1)
            if ~isempty(C{i,1})
                Cluster{i,1}(1,:) = C{i,1}(1,:); Cluster{i,1}(2,:) = C{i,1}(2,:);
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
            rep = {CNPSummary(:).Sample};
            ORI = {CNPSummary(:).ORI_Crossing};
        end
        
    elseif IndataType == "Adist"
        for i=1:size(Adist,1)
            if ~isempty(Adist{i,1})
                Cluster{i,1}(1,:) = Adist{i,1}(:,2)'; Cluster{i,1}(2,:) = Adist{i,1}(:,3)';
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
            rep = {CNPSummary(:).Sample};
            ORI = {CNPSummary(:).ORI_Crossing};
        end
        
    elseif IndataType == "denovo"     
        for i=1:size(denovo,1)
            if ~isempty(denovo{i,1})
                Cluster{i,1}(1,:) = denovo{i,1}(1,:)'; Cluster{i,1}(2,:) = denovo{i,1}(1,:)';
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
            rep = {CNPSummary(:).Sample};
            ORI = {CNPSummary(:).ORI_Crossing};
        end
    end
    
elseif ismember(IndataType,["del","dup"])
    D=load(IndataPath + Indata);
    DD.del=D.deldup.del(dataset);
    DD.dup=D.deldup.dup(dataset);
    delstart={DD.del(:).start}';
    deledge={DD.del(:).edge}';
    
    if IndataType == "del"
        for i=1:size(DD.del,2)
            if ~isempty(delstart{i,1})
                Cluster{i,1}(1,:) = delstart{i,1}(:)'; Cluster{i,1}(2,:) = deledge{i,1}(:)';
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
            rep = {DD.del(:).sample};
        end
        
    elseif IndataType == "dup"       
        for i=1:size(DD.dup,2)
            if ~isempty(dupstart{i,1})
                Cluster{i,1}(1,:) = dupstart{i,1}(:)'; Cluster{i,1}(2,:) = dupedge{i,1}(:)';
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
            rep = {DD.dup(:).sample};
        end
    end
    
elseif IndataType == "txtListStart"
    fid = fopen(IndataPath + Indata);
    imp = textscan(fid,'%f %s %s','delimiter','\t');
    fclose(fid);
    txt.S=imp{1}; txt.E=imp{1}; txt.Ref = imp{2}; txt.Var = imp{3}; 
    clear imp
 
        for i=1:numel(txt.S)
            Cluster{i,1}(1,:) = txt.S(i); Cluster{i,1}(2,:) = txt.E(i);
            rep{i} = string(1);
            mut_ref{i} = txt.Ref(i);
            mut_var{i} = txt.Var(i);
        end
        
elseif IndataType == "txtListStartEnd"
    fid = fopen(IndataPath + Indata);
    imp = textscan(fid,'%s %s %f %s %s %s %s  %s %s  %s %s %s','delimiter','\t');
    fclose(fid);
    txt.S=imp{9}; txt.E=imp{10}; %txt.Ref = imp{2}; txt.Var = imp{3}; 
    clear imp
 
        for i=1:numel(txt.S)
            Cluster{i,1}(1,:) = txt.S(i); Cluster{i,1}(2,:) = txt.E(i);
            rep{i} = string(1);
            %mut_ref{i} = txt.Ref(i);
            %mut_var{i} = txt.Var(i);
        end
end


% in the case that we are looking at denovos, deletions and duplications, acc. genes
% are not excluded!
if ismember(IndataType,["del","dup","denovo"])
        exclude_thr = 1.1;
end
%% HitOutput is generated
%Here for every replicate i all cluster c are checked for hits with genes and all different hits are written as entry into HitOutput.
%   -> This means: genes can appear multiple times per replicate -> this is
%      fixed in GeneOutput
HitOutput = struct('RepNo',[],'Cluster',[],'Type',[],'Frac',[],'locus_tag',[],'Name',[],'product',num2cell(repmat("",1,1000)),'GeneEndePos',[],'GeneDirection',[],'MutationPos',[],...
    'Ref',[],'Var',[],'Sample',[]);
m = 0;
for i=1:size(Cluster,1) % loops over all samples
    % fix the problem with ORI first
    d = 1;
    while d <= size(Cluster{i},2) % loops over all cluster in ith sample
        % if we have a cluster where start > edge, we assume ORI crossing and
        % split up this entry in 2: start(1):recipsize & 1:end(1)
        if Cluster{i}(1,d) > Cluster{i}(2,d)
            first = [1 Cluster{i}(2,1)]; last = [Cluster{i}(1,1) recipsize];
            Cluster{i}(:,1) = first; Cluster{i}(:,end+1) = last;
            disp('A cluster was found that crossed the ORI. It is cut in 2 pieces!');
        end
        d = d+1;
    end
    
    for c=1:size(Cluster{i},2) % loops over all cluster in ith sample
        ptail    = find(((Cluster{i}(1,c) > [bedFile_tmp.start]) & (Cluster{i}(1,c) < [bedFile_tmp.ende]) & (Cluster{i}(2,c) >= [bedFile_tmp.ende])));
        complete = find((Cluster{i}(1,c) <= [bedFile.start]) & (Cluster{i}(2,c) >= [bedFile.ende]));
        pin      = find((Cluster{i}(1,c) > [bedFile.start]) & (Cluster{i}(2,c) < [bedFile.ende]));
        pfront   = find((Cluster{i}(1,c) <= [bedFile.start]) & (Cluster{i}(2,c) > [bedFile.start]) & (Cluster{i}(2,c) < [bedFile.ende]));
       
        allIdx   = [ptail complete pin pfront];
        allTypes = [repmat("tail",1,numel(ptail)) repmat("complete",1,numel(complete))...
            repmat("in",1,numel(pin)) repmat("front",1,numel(pfront))];
        
        for t=1:numel(allIdx) 
            idxH = allIdx(t);
            m = m+1;            
            HitOutput(m).Type     = allTypes(t);
            HitOutput(m).RepNo    = i; 
            HitOutput(m).Sample   = rep{i};    
            HitOutput(m).Cluster  = c;
            HitOutput(m).Name     = string(bedFile(idxH).Name);
            HitOutput(m).locus_tag= string(bedFile(idxH).locus_tag); 
            if ~isempty(bedFile(idxH).product)
            HitOutput(m).product  = string(bedFile(idxH).product); 
            end
            HitOutput(m).MutationPos = Cluster{i}(1,c);
            HitOutput(m).Ref = string(mut_ref{i});
            HitOutput(m).Var = string(mut_var{i});
            HitOutput(m).GeneEndePos  = bedFile(idxH).ende;
            HitOutput(m).GeneDirection = string(bedFile(idxH).direction);
            
            if allTypes(t) == "tail"                                                
                HitOutput(m).Frac = (bedFile(idxH).ende - Cluster{i}(1,c) + 1)/bedFile(idxH).L;
            end

            if allTypes(t) == "complete"                
                HitOutput(m).Frac = 1;
            end

             if allTypes == "in"                             
                HitOutput(m).Frac = (Cluster{i}(2,c) - Cluster{i}(1,c) + 1)/bedFile(idxH).L;
            end

            if allTypes == "front"                              
                HitOutput(m).Frac = (Cluster{i}(2,c) - bedFile(idxH).start + 1)/bedFile(idxH).L;
            end
        end
        
        clear complete pin pfront ptail c co d pf pi pt
    end
end

HitOutput = HitOutput(1:numel([HitOutput.RepNo]))
%% GeneOutput and RepSummary are generated
GeneOutput = struct('RepNo',[],'Cluster',[],'Type',[],'Frac',[],'locus_tag',[],'Name',[],'product',[],'GeneEndePos',[],'GeneDirection',[],'MutationPos',[],...
    'Ref',[],'Var',[],'Sample',[]);
RepSummary = struct('RepNo',[],'locus_tag_Hit',[],'NameHit',[],'Sample',[]);

l = 0;
for i=1:numel(rep)
    idx = find([HitOutput.RepNo] == i);
    LTcollectRep = [HitOutput(idx).locus_tag];
    for j=1:numel(idx)
        match = find(strcmp(LTcollectRep(j),LTcollectRep));
        matchD = match(match~=j);
        if length(match)==1
            l = l + 1;
            GeneOutput(l) = HitOutput(idx(j));
        elseif matchD>j
            l = l + 1;
            GeneOutput(l) = HitOutput(idx(j));
            %add the other gene infos to this entry
            GeneOutput(l).Frac = sum([HitOutput(idx(match)).Frac]);
            GeneOutput(l).Cluster = [HitOutput(idx(match)).Cluster];
            GeneOutput(l).Type = {HitOutput(idx(match)).Type};
        end
    end
    
    %Now for each Replicate I collect the summary
    idxUniq = find([GeneOutput.RepNo] == i);
    RepSummary(i).locus_tag_Hit  = [GeneOutput(idxUniq).locus_tag];  
    RepSummary(i).NameHit        = [GeneOutput(idxUniq).Name];     
    RepSummary(i).RepNo = i;    
    RepSummary(i).Sample = rep{i};
    RepSummary(i).GeneHitNo = length(idxUniq);
end

% %% Multiple HitStatistics: How often is each Gene hit over all replicates??
% % Prepare gene list so that it excludes the genes that are excAccmm and
% %that are above the threshold
% acm = {{excAccmm(:).BSU}'}; accmm.BSU = vertcat(acm{:});
% 
% mask = [excAccmm(:).FracMean]>=exclude_thr;
% ex = {{excAccmm(mask).BSU}'}; exc.BSU = vertcat(ex{:});
% 
% mlSNPsinGene = 0;
% MultiHitStat = struct('BSU',[],'Genename',[],'SumHit',[],'HitIndex',[],'FracMean',[],'mlGeneIdent',[],'Samples',[],'accxluded',[],'FracList',[]);
% c = {GeneOutput(:).BSU}';          BSUcollect = vertcat(c{:});
% Samplecollect = {GeneOutput(:).Sample}';
% n = 0;
% % here, i goes through all the entries in the bed file
% for i=1:numel(bedFile_tmp.BSU)
%     
%     MultiHitStat(i).BSU = bedFile_tmp.BSU{i};
%     MultiHitStat(i).Genename = bedFile_tmp.GN{i};
%     
%     %add the gene Identity of the i-th gene based on the masterlist to
%     %MultiHitStat- if the i-th gene also appears on the excAnnMM list, then
%     %calculate it a bit differently
%     idx_acc = find(strcmp(bedFile_tmp.BSU{i},accmm.BSU));
%     if ~isempty(idx_acc)
%         % here a very small number is added just to not to devide by 0
%         l = bedFile_tmp.L(i)*(1-(excAccmm(idx_acc).FracMean));
%         snps = sum(ismember(refchr.pos,bedFile_tmp.S(i):bedFile_tmp.E(i)));
%         MultiHitStat(i).mlGeneIdent =  1-(snps/l);
%         if snps>l; snps = l; MultiHitStat(i).mlGeneIdent =  1-(snps/l); end
%         if l == 0; l=0.0000001; MultiHitStat(i).mlGeneIdent =  0; end
%     else
%         MultiHitStat(i).mlGeneIdent =  1-(sum(ismember(refchr.pos,bedFile_tmp.S(i):bedFile_tmp.E(i)))/bedFile_tmp.L(i));
%     end
%     
%     % if the gene also appears in exc.BSU, mark it as excluded in the
%     % accxluded field
%     exc_idx = find(strcmp(bedFile_tmp.BSU{i},exc.BSU));
%     % at which index idx did the i-th gene fit the geneoutput list?
%     idx = find(strcmp(bedFile_tmp.BSU{i},BSUcollect));
%     if ~isempty(exc_idx)
%         MultiHitStat(i).accxluded = 1;
%         MultiHitStat(i).SumHit = 0;
%         MultiHitStat(i).FracMean = 0;
%         
%         % if the gene was hit at least once, calculate ..
%     elseif isempty(exc_idx) && ~isempty(idx)
%         MultiHitStat(i).accxluded = 0;
%         MultiHitStat(i).SumHit = numel(idx);
%         MultiHitStat(i).HitIndex = idx;
%         MultiHitStat(i).FracList = [GeneOutput(idx).Frac];
%         MultiHitStat(i).FracMean = mean([GeneOutput(idx).Frac]);
%         MultiHitStat(i).Samples = {Samplecollect{idx}};
%         
%     elseif isempty(exc_idx) && isempty(idx)
%         MultiHitStat(i).accxluded = 0;
%         MultiHitStat(i).SumHit = 0;
%         MultiHitStat(i).FracMean = 0;
%     end
%     
%     clear idx
% end
% 
% %% Create a list of the most interesting hotspots
% % save data
% saveit = input('Do you want to save the HotColdGenes data? (yes/no) ');
% if isempty(saveit) || saveit == 1
%     cutit = input('Which minimal number of MultiHits interests you? Set cutoff: ');
%     cutoff=cutit; c = 0; Liste = struct('BSU',[],'Genename',[],'SumHit',[],'HitIndex',[],'FracMean',[],'mlGeneIdent',[],'Samples',[],'accxluded',[],'FracList',[]);
%     for i=1:numel(bedFile_tmp.BSU)
%         if MultiHitStat(i).SumHit>=cutoff
%             c = c+1;
%             Liste(c) = MultiHitStat(i);
%         end
%     end
%     
%     if strcmp(which,'denovo')
%         Liste = rmfield(Liste,{'HitIndex','mlGeneIdent','FracMean','accxluded','FracList'});
%     else
%         Liste = rmfield(Liste,{'HitIndex','mlGeneIdent'});
%     end
%     
%     T=struct2table(Liste);
%     writetable(T,[savepath 'den2Genes_Wscy20.txt'],'Delimiter',' ')
% else
%     disp('Ok, HotColdGenes are not saved.')
% end
% 
% %% Hot and cold plot
% % be aware that the genes that appear here as acc. with ident=0 are the ones that
% % are excluded with the threshold
% 
% rows = 6; totalgenes = numel(bedFile_tmp.L);
% ident = [MultiHitStat(:).mlGeneIdent]; sumhit = [MultiHitStat(:).SumHit];
% left_color = [0.1 0.1 0.1];
% right_color = 'b';
% 
% figure(1)
% subplot(6,1,1)
% title('Hot and Cold plot')
% hold on
% for i=1:rows
%     subplot(rows,1,i)
%     region = [floor(1+(i-1)*(totalgenes/rows)):ceil((totalgenes/rows)*i)];
%     yyaxis left;
%     plot(region,ident(region),'Color',[left_color 0.3])
%     ylim([0.8 1])
%     hold on
%     yyaxis right
%     bar(region,sumhit(region),'BarWidth',1,'FaceColor',right_color)
%     xlim([1+(i-1)*(totalgenes/rows) (totalgenes/rows)*i ])
%     ylim([0 max(sumhit)+1])
%     ax = gca;
%     ax.YAxis(1).Color = left_color;
%     ax.YAxis(2).Color = right_color;
%     
% end
% subplot(6,1,ceil(rows/2))
% yyaxis left; ylabel('Identity')
% yyaxis right; ylabel('# of hits')
% subplot(6,1,rows)
% xlabel('Gene number');
% 
% %% MultiHitStatistic Plot - raw data
% if strcmp(which,'txt') ~= 1
%     % in this plot the acc. genes that you defined with exclude_thr are
%     % excluded
%     sumhitmask_woacc = [MultiHitStat(:).accxluded]==0;
%     sumhit_woacc = sumhit(sumhitmask_woacc);
%     
%     figure(2);
%     h1=histogram(sumhit,'BinWidth',1); hold on
%     h1.BinEdges = [h1.BinEdges max(h1.BinEdges)+1] - h1.BinWidth/2;
%     h2=histogram(sumhit_woacc,'BinWidth',1,'EdgeColor','r','FaceColor','none');
%     h2.BinEdges = [h2.BinEdges max(h2.BinEdges)+1] - h2.BinWidth/2;
%     set(gca, 'YScale', 'log')
%     title('Multiple Hit Statistic')
%     xlabel('Gene Hits'); ylabel('# of hits')
%     xlim([-1 numel(h1.Values)+2]); ylim([1 max(h1.Values)+500])
%     legend([h1 h2],{'Ws cy20','Ws cy20 wo AccGenes'})
% end
% 
% allBSU = {bedFile_tmp.BSU}'; allBSUU = vertcat(allBSU{:});
% for i=1:numel([MultiHitStat(:).FracMean])
%     idx=find(strcmp(string(MultiHitStat(i).BSU),allBSUU));
%     MultiHitStat(i).GeneNum = idx;
% end


save("/home/isabel/Dokumente/P5_ExpEvol_Ngo/DNA/analysis/Mutations_2_Genes/" + "Ngo2MS11_Mutations_2_GeneOutput_wIndels.mat",'GeneOutput')
% as table in .txt
T_CDS=struct2table(GeneOutput);
writetable(T_CDS,"/home/isabel/Dokumente/P5_ExpEvol_Ngo/DNA/analysis/Mutations_2_Genes/" + "Ngo2MS11_Mutations_2_GeneOutput_wIndels.txt",'Delimiter','\t')