%
%
%%%%%%%%%%%%%%%%   CNP2genes created by MonaIsa  %%%%%%%%%%%%%%%%%%%%%%%%
%
%  What this script does: 
%     -- detects for all genomic replacements, detected as clusters, which genes were affected by this.
%     -- counts how often genes are hit and does Multi-Hit Statistics for them
%
%%%%%%%%%%%%%%
%  Input: - Cluster information: either CNPSummary, deldup or a txt file
%         - bed file; masterlist
%         - dataset: where in SNPSummary or deldup are the data of
%         interest?
%         - savepath: where are the data saved?
%         - for the Multi-Hit Statistics the accessory and mutimapper genes are excluded, therefore we need
%           - Input: AccMM2Genes_Hits.mat
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

clear all; 
close all

%% Here, you have to make some inputs:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Data Input   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% different data inputs, choose one:
IndataType = "Adist";       % from CNPSummary/ CNPISummary: "Adist" or "denovo" or "SPI"
                            % from CNPISummary: "indel"
                            % from Output_DelDup: "dup" or "del"
                            % from indelSummary: "indel"
                            % plain txt file with either start and end: "txtListStartEnd" (start end length [%f %f %f])
                            % or only 1 position (SNP): "txtListSNP"
% Here give the data file
% either: CNPSummary.mat, Output_deldup.mat or txt document

IndataPath = "../IN";
Indata     = "CNPSummary_Wns17.mat";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Additional Files  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

donor = "Bval"; % "Bspiz", "Bval", "Batro"

pathLists = "/DetectionGV/dictionaries_Bacillus/";
if strcmp(donor, "Bval")
    % Bval donor
    masterlist = pathLists + "ml/" + "mlBval2Bs166NCe_v1.txt";
    accgenome  = pathLists + "acc/" + "Bval_AccMM2Genes.mat";
elseif strcmp(donor, "Bspiz")
    % Bspiz donor
    masterlist = pathLists + "ml/" + "mlW232Bs166NCe_v1.txt";
    accgenome  = pathLists + "acc/" +"W23_AccMM2Genes.mat";
elseif strcmp(donor, "Batro")
    % BA donor
    masterlist = pathLists + "ml/" + "mlBatro2Bs166NCe_v1.txt";
    accgenome  = pathLists + "acc/" +"Batro_AccMM2Genes.mat";
elseif strcmp(donor, "Geo")
    % Geo donor
    masterlist = pathLists + "ml/" + "mlGeo2Bs166NCe_v1.txt";
    accgenome  = pathLists + "acc/" +"Batro_AccMM2Genes.mat"; %
else
    error("Something went wrong with your donor declaration.. Is your donor really %s?", donor);
end


% now we need the information, where which gene is in your organism
recipbed   = "/DetectionGV/dictionaries_Bacillus/Bs166NCe_June2021.bed.mat";

% Do you want to only search for subset of genes given in the bed?
SearchInSpecGenes = "OFF";
specGenes         = "....txt";

%%%%%%%     Load variables   %%%%%%%%%%%%
% which accmm genes to exclude? exclude_thr sets to which frac a gene must
% be accmm to be excluded; exclude_thr sets the lower limit
% % % if exclude_thr = 1, then only genes, that are complete accmm are exc.
% % % if exclude_thr = 1.1, then no gene is excluded
exclude_thr = 1.1;

recipsize = 4215607;

% the plots and hotcolds are saved here
savepath = "../OUT";

% Do you want to save the HotColdGenes data?
saveHotColdData = "OFF";
% Which minimal number of MultiHits interests you? Set cutoff:
cutoff = 1;

%% Here the data is loaded

%Load recipient/donor specific master list
fid = fopen(masterlist); imp = textscan(fid,'%f %s %s');
fclose(fid);
refchr.pos=imp{1}; refchr.atcg=imp{3};
clear imp
clear prompt str fid ans imp

% % load the acc and mm genes for MultiHitStat
excAccmm = load(accgenome);
excAccmm = excAccmm.AccMM2Genes;

% %Load recipient annonated bed file
bed     = load(recipbed);
bedFile_tmp = bed.uniqList;

% in case you are only interested in the occurrence of specific genes ..
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
    bedFile = bedFile_tmp(speIdx);

else % take the complete bed file
    bedFile = bedFile_tmp;
end

% add the gene length information to bed
tmp = struct("L", num2cell([bedFile.ende] - [bedFile.start] +1));
[bedFile.L] = deal(tmp.L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read in main data of interest %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(IndataType,["Cdist","Adist","denovo","indel"])
    CNPSummary=load(IndataPath + Indata);
    if isempty(dataset)
        CNPSummary=CNPSummary.CNPSummary;
    else
        data_mask = contains([CNPSummary.CNPSummary.Sample],dataset);
        CNPSummary=CNPSummary.CNPSummary(data_mask);
    end
    % this is same for all
    C      = {CNPSummary(:).C}';
    Adist  = {CNPSummary(:).Adist}';
    denovo = {CNPSummary(:).denovo}';
    indel  = {CNPSummary(:).indel}';
    rep    = {CNPSummary(:).Sample};
    ORI    = {CNPSummary(:).ORI_Crossing};

    if IndataType == "Cdist"
        for i=1:size(C,1)
            if ~isempty(C{i,1})
                Cluster{i,1}(1,:) = C{i,1}(1,:); Cluster{i,1}(2,:) = C{i,1}(2,:);
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end

        end

    elseif IndataType == "Adist"
        for i=1:size(Adist,1)
            if ~isempty(Adist{i,1})
                Cluster{i,1}(1,:) = Adist{i,1}(:,2)'; Cluster{i,1}(2,:) = Adist{i,1}(:,3)';
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
        end

    elseif IndataType == "denovo"
        for i=1:size(denovo,1)
            if ~isempty(denovo{i,1})
                de_notClu = denovo{i,1}.pos(2,:)' ~= 1;
                Cluster{i,1}(1,:) = denovo{i,1}.pos(1,de_notClu)'; Cluster{i,1}(2,:) = denovo{i,1}.pos(1,de_notClu)';
                Cluster{i,2}(1,:) = denovo{i,1}.ref(de_notClu)'; Cluster{i,3}(1,:) = denovo{i,1}.var(de_notClu)';
               
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
        end

    elseif IndataType == "indel"
        for i=1:size(indel,1)
            if ~isempty(indel{i,1})
                ind_notClu = denovo{i,1}.pos(2,:)' ~= 1;
                Cluster{i,1}(1,:) = indel{i,1}.pos(1,ind_notClu)'; Cluster{i,1}(2,:) = indel{i,1}.pos(1,ind_notClu)';
                Cluster{i,2}(1,:) = indel{i,1}.ref(ind_notClu)'; Cluster{i,3}(1,:) = indel{i,1}.var(ind_notClu)';
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
        end

    elseif IndataType == "SPI"
        for i=1:size(SPI,1)
            if ~isempty(SPI{i,1})
                Cluster{i,1}(1,:) = SPI{i,1}(:,1); Cluster{i,1}(2,:) = SPI{i,1}(:,1);
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
        end
    end


elseif ismember(IndataType,["del","dup"])
    D=load(IndataPath + Indata);
    data_mask = contains([D.deldup.del.sample],dataset);
    if isempty(dataset)
        DD.del=D.deldup.del;
        DD.dup=D.deldup.dup;
    else
        DD.del=D.deldup.del(data_mask);
        DD.dup=D.deldup.dup(data_mask);
    end

    delstart={DD.del(:).start}';
    deledge={DD.del(:).edge}';
    dupstart={DD.dup(:).start}';
    dupedge={DD.dup(:).edge}';

    if IndataType == "del"
        for i=1:size(DD.del,2)
            if ~isempty(delstart{i,1})
                Cluster{i,1}(1,:) = delstart{i,1}(:)'; Cluster{i,1}(2,:) = deledge{i,1}(:)';
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
        end
        rep = {DD.del(:).sample};

    elseif IndataType == "dup"
        for i=1:size(DD.dup,2)
            if ~isempty(dupstart{i,1})
                Cluster{i,1}(1,:) = dupstart{i,1}(:)'; Cluster{i,1}(2,:) = dupedge{i,1}(:)';
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
        end
        rep = {DD.dup(:).sample};

    end

elseif IndataType == "txtListStartEnd"
    fid = fopen(IndataPath + Indata);
    imp = textscan(fid,'%f %f %f','delimiter',' ');
    fclose(fid);
    txt.S=imp{1}; txt.E=imp{2};
    clear imp

    for i=1:numel(txt.S)
        Cluster{i,1}(1,:) = txt.S(i); Cluster{i,1}(2,:) = txt.E(i);
        rep{i} = string(1);

    end
end


% when looking at denovos, del, dup -> acc. genes not excluded!
if ismember(IndataType,["del","dup","denovo", "indel"])
    exclude_thr = 1.1;
end

sampleNum = size(Cluster,1);

%% HitOutput is generated

%Here for every replicate i all cluster c are checked for hits with genes and all different hits are written as entry into HitOutput.
%   -> This means: genes can appear multiple times per replicate -> this is
%      fixed in GeneOutput

HitOutput = struct('RepNo',[],'Cluster',[],'Type',[],'Frac',[],'locustag',[],'Name',[],'mlHitIdent',[], 'product',[],'GeneEndePos',[],...
    'GeneDirection',[],'Sample',[],'ChangeStart',[],'ChangeEnd',[],'ChangeRef',[],'ChangeVar',[]);
m = 0;
for i=1:sampleNum % loops over all samples
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
        ptail    = find(((Cluster{i}(1,c) > [bedFile.start]) & (Cluster{i}(1,c) < [bedFile.ende]) & (Cluster{i}(2,c) >= [bedFile.ende])));
        complete = find((Cluster{i}(1,c) <= [bedFile.start]) & (Cluster{i}(2,c) >= [bedFile.ende]));
        pin      = find((Cluster{i}(1,c) > [bedFile.start]) & (Cluster{i}(2,c) < [bedFile.ende]));
        pfront   = find((Cluster{i}(1,c) <= [bedFile.start]) & (Cluster{i}(2,c) > [bedFile.start]) & (Cluster{i}(2,c) < [bedFile.ende]));

        allIdx   = [ptail complete pin pfront];
        allTypes = [repmat("tail",1,numel(ptail)) repmat("complete",1,numel(complete))...
            repmat("in",1,numel(pin)) repmat("front",1,numel(pfront))];

        for t=1:numel(allIdx)
            idxH = allIdx(t);
            m = m+1;
            if any(strcmp(IndataType, ["denovo","indel"]))
                HitOutput(m).ChangeRef  = Cluster{i,2}(c);
                HitOutput(m).ChangeVar  = Cluster{i,3}(c);
            end
            HitOutput(m).ChangeStart= Cluster{i}(1,c);
            HitOutput(m).ChangeEnd  = Cluster{i}(2,c);
            
            HitOutput(m).Type       = allTypes(t);
            HitOutput(m).RepNo      = i;
            HitOutput(m).Sample     = rep{i};
            HitOutput(m).Cluster    = c;
            HitOutput(m).Name       = string(bedFile(idxH).Name);
            HitOutput(m).locustag   = string(bedFile(idxH).locus_tag);
            HitOutput(m).product    = string(bedFile(idxH).product);
            HitOutput(m).GeneEndePos = bedFile(idxH).ende;
            HitOutput(m).GeneDirection = string(bedFile(idxH).direction);

            if allTypes(t) == "tail"
                HitOutput(m).Frac = (bedFile(idxH).ende - Cluster{i}(1,c) + 1)/bedFile(idxH).L; % Cluster{i}(1,c) <- start cluster
                snps        = sum(ismember(refchr.pos,Cluster{i}(1,c):bedFile(idxH).ende)); % Cluster{i}(2,c) <- end cluster
            end

            if allTypes(t) == "complete"
                HitOutput(m).Frac = 1;
                snps        = sum(ismember(refchr.pos,bedFile(idxH).start:bedFile(idxH).ende));
            end

            if allTypes(t) == "in"
                HitOutput(m).Frac = (Cluster{i}(2,c) - Cluster{i}(1,c) + 1)/bedFile(idxH).L;
                snps        = sum(ismember(refchr.pos,Cluster{i}(1,c):Cluster{i}(2,c)));
            end

            if allTypes(t) == "front"
                HitOutput(m).Frac = (Cluster{i}(2,c) - bedFile(idxH).start + 1)/bedFile(idxH).L;
                snps        = sum(ismember(refchr.pos,bedFile(idxH).start:Cluster{i}(2,c)));
            end

            % Calculate the identity of the hit part of the gene
            % (i) count snps (in if loops - depends on hitType)
            % (ii) divide by hit length (fracHit * total length)
            % (iii) subtract from 1 to have ident not divergence
            HitOutput(m).mlHitIdent = 1 - snps/(HitOutput(m).Frac*bedFile(idxH).L);

        end

        clear complete pin pfront ptail c co d pf pi pt snps
    end
end

%% GeneOutput and RepSummary are generated

GeneOutput = struct('RepNo',[],'Cluster',[],'Type',[],'Frac',[],'locustag',[],'Name',[],'mlHitIdent',[], 'product',[],'GeneEndePos',[],...
    'GeneDirection',[],'Sample',[],'ChangeStart',[],'ChangeEnd',[],'ChangeRef',[],'ChangeVar',[]);
RepSummary = struct('RepNo',[],'locustag_Hit',[],'NameHit',[],'Sample',[]);

l = 0;
for i=1:numel(rep)
    idx = find([HitOutput.RepNo] == i);
    LTcollectRep = [HitOutput(idx).locustag];
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
            GeneOutput(l).mlHitIdent = [HitOutput(idx(match)).mlHitIdent];
        end
    end

    %Now for each Replicate, collect the summary
    idxUniq = find([GeneOutput.RepNo] == i);
    RepSummary(i).locustag_Hit  = [GeneOutput(idxUniq).locustag];
    RepSummary(i).NameHit        = [GeneOutput(idxUniq).Name];
    RepSummary(i).RepNo = i;
    RepSummary(i).Sample = rep{i};
    RepSummary(i).GeneHitNo = length(idxUniq);
end

%% Multiple HitStatistics: How often is each Gene hit over all replicates??

mask = [excAccmm.FracMean]>=exclude_thr;
exc  = [excAccmm(mask).locustag];

% mlSNPsinGene = 0;
MultiHitStat = struct('locustag',[],'GN',[],'SumHit',[],'HitIndex',[],'FracMean',[],'mlGeneIdent',[],'Samples',[],'accxluded',[],'FracList',[]);

% here, i goes through all the entries in the bed file
for i=1:numel([bedFile.locus_tag])

    MultiHitStat(i).locustag = bedFile(i).locus_tag;
    MultiHitStat(i).GN       = bedFile(i).Name;

    % add the gene Identity of the i-th gene based on the masterlist to
    % MultiHitStat
    % % % --- if the i-th gene also appears on the excAnnMM list, then
    % % %     calculate it a bit differently

    idx_acc = find(strcmp(bedFile(i).locus_tag,[excAccmm.locustag]));
    if ~isempty(idx_acc)

        % here a very small number is added just to not to devide by 0
        l_bp_Acc    = bedFile(i).L * excAccmm(idx_acc).FracMean;
        l_bp_nonAcc = bedFile(i).L - l_bp_Acc;
        snps        = sum(ismember(refchr.pos,bedFile(i).start:bedFile(i).ende));

        MultiHitStat(i).mlGeneIdent     =  (1 - excAccmm(idx_acc).FracMean) * (1-(snps/l_bp_nonAcc)); % Ident of complete gene = FracNonAcc * Ident_NonAcc (+ FracAcc * Ident_Acc ) <- Ident_Acc = 0 !


        if snps > l_bp_nonAcc
            snps = l_bp_nonAcc;
            MultiHitStat(i).mlGeneIdent =  1-(snps/l_bp_nonAcc);
        end

        if l_bp_nonAcc == 0
            MultiHitStat(i).mlGeneIdent =  0;
        end

    else
        MultiHitStat(i).mlGeneIdent =  1-(sum(ismember(refchr.pos,bedFile(i).start:bedFile(i).ende))/bedFile(i).L);
    end

    % if the gene also appears in exc.BSU, mark it as excluded in the
    % accxluded field
    exc_decide = find(strcmp(bedFile(i).locus_tag,exc));
    % at which index idx did the i-th gene fit the geneoutput list?
    idx = find(strcmp(bedFile(i).locus_tag,[GeneOutput.locustag]));

    if ~isempty(exc_decide)
        MultiHitStat(i).accxluded = 1;
        MultiHitStat(i).SumHit = 0;
        MultiHitStat(i).FracMean = 0;

        % if the gene was hit at least once, calculate ..
    elseif isempty(exc_decide) && ~isempty(idx)

        MultiHitStat(i).accxluded = 0;
        MultiHitStat(i).SumHit    = numel(idx);
        MultiHitStat(i).HitIndex  = idx;
        MultiHitStat(i).FracList  = [GeneOutput(idx).Frac];
        MultiHitStat(i).FracMean  = mean([GeneOutput(idx).Frac]);
        MultiHitStat(i).Samples   = [GeneOutput(idx).Sample];

    elseif isempty(exc_decide) && isempty(idx)
        MultiHitStat(i).accxluded = 0;
        MultiHitStat(i).SumHit    = 0;
        MultiHitStat(i).FracMean  = 0;
    end

    clear idx
end

%% Create a list of the most interesting hotspots
% save data
if saveHotColdData == "ON"
    c = 0;
    Liste = struct('locustag',[],'Genename',[],'SumHit',[],'HitIndex',[],'FracMean',[],'mlGeneIdent',[],'Samples',[],'accxluded',[],'FracList',[]);
    for i=1:numel([bedFile.locus_tag])
        if MultiHitStat(i).SumHit>=cutoff
            c = c+1;
            Liste(c) = MultiHitStat(i);
        end
    end

    if strcmp(which,'denovo')
        Liste = rmfield(Liste,{'HitIndex','mlGeneIdent','FracMean','accxluded','FracList'});
    else
        Liste = rmfield(Liste,{'HitIndex','mlGeneIdent'});
    end

    T=struct2table(Liste);
    writetable(T,[savepath 'den2Genes_Wscy20.txt'],'Delimiter',' ')
else
    disp('Ok, HotColdGenes are not saved.')
end

%% Here, some plots for hot spots and multihit statistics

% % % Hot and cold plot
% be aware that the genes that appear here as acc. with ident=0 are the ones that
% are excluded with the threshold

rows = 6; totalgenes = numel([bedFile.L]);
ident = [MultiHitStat(:).mlGeneIdent];
sumhit = [MultiHitStat(:).SumHit];
left_color = [0.1 0.1 0.1];
right_color = 'b';

figure(1)
subplot(6,1,1)
title('Hot and Cold plot')
hold on
for i=1:rows
    subplot(rows,1,i)
    region = [floor(1+(i-1)*(totalgenes/rows)):ceil((totalgenes/rows)*i)];
    yyaxis left;
    plot(region,ident(region),'Color',[left_color 0.3])
    ylim([0.8 1])
    hold on
    yyaxis right
    bar(region,sumhit(region),'BarWidth',1,'FaceColor',right_color)
    xlim([1+(i-1)*(totalgenes/rows) (totalgenes/rows)*i ])
    ylim([0 max(sumhit)+1])
    ax = gca;
    ax.YAxis(1).Color = left_color;
    ax.YAxis(2).Color = right_color;

end
subplot(6,1,ceil(rows/2))
yyaxis left; ylabel('Identity')
yyaxis right; ylabel('# of hits')
subplot(6,1,rows)
xlabel('Gene number');

% % % MultiHitStatistic Plot - raw data

% in this plot the acc. genes that you defined with exclude_thr are
% excluded
sumhitmask_woacc = [MultiHitStat(:).accxluded]==0;
sumhit_woacc     = sumhit(sumhitmask_woacc);

figure(2);
h1=histogram(sumhit,'BinWidth',1); hold on
h1.BinEdges = [h1.BinEdges max(h1.BinEdges)+1] - h1.BinWidth/2;
h2=histogram(sumhit_woacc,'BinWidth',1,'EdgeColor','r','FaceColor','none');
h2.BinEdges = [h2.BinEdges max(h2.BinEdges)+1] - h2.BinWidth/2;
set(gca, 'YScale', 'log')
title('Multiple Hit Statistic')
xlabel('Gene Hits'); ylabel('# of hits')
xlim([-1 numel(h1.Values)+2]); ylim([1 max(h1.Values)+500])
legend([h1 h2],{'cy20','cy20 wo AccGenes'})


