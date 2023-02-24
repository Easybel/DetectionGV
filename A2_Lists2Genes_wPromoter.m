%% the idea of this script is to find promoter regions
%
%
% in/ out:
% % %  -- Cluster data (or denovo, deldup, list of regions ...)
% % %  -- bed file of the species in question
% % %  -- List of promoters of Bsub
% % %  -- masterlist

% % what the script does:
% --- this script searches for all genes that were hit by HGT (clusters, denovos, deldup ...)
% --> out: HitOutput (here, genes can appear more then once if they were hit by different clusters)
% --- then it goes through HitOutput and brings all hits together per gene
% --> out: GeneOutput
%
% ----- Additionally: this script can, if PROM == "ON", go through all promoters
% ----- and check if they are hit
% ----> out: HitOutputProm
% ----- then we check for each gene in GeneOutput if the according promoter
%       was hit
% promoter was hit
% whether the promoter of these genes were also hit


close all; clear all
%% set paths and read in data
% basePath = "/home/isabel/sciebo/ResultsShared/kleinesPaper/";
basePath    = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/";
operonPath  = basePath + "allLists/" + "Bs166NC_OperonsSubtiwiki_20210903.mat";
recipbed   = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/Bs166NCe_June2021.bed.mat";
pseudoGenes = basePath + "allLists/" + "Bs166NC_PseudoGenesSubtiwiki_20210902.mat";
% Here give the data file
% either: CNPSummary.mat, Output_deldup.mat or txt document
%IndataPath  = "/home/isabel/sciebo/ResultsShared/DFE_HighThroughput/4_DNASeq/";
IndataPath  = "/home/isabel/Documents/sciebo/DFE_fromSciebo/2022_withSelection/DNAseq/";


% Do you want to only search for subset of genes given in the bed?
SearchInSpecGenes = "OFF";
specGenes         = "/home/isabel/Documents/Doktorarbeit_Mai2022/P1_EvolExp_CLASSIC_Bacillus/0_dictionaries/NCIB3610_wildstrain/Mapped_Bs166_2_NCIB3610/Bs166_2_NCIB3610_HitOutput_NonSnyBSU.txt";

% different data inputs, choose one:
IndataType = "Adist";       % from CNPSummary/ CNPISummary: "Adist" or "denovo" or "SPI"
                            % from CNPISummary: "indel"
                            % from Output_DelDup: "dup" or "del"
                            % from indelSummary: "indel"
                            % plain txt file with either start and end: "txtListStartEnd" (start end length [%f %f %f])
                            % or only 1 position (SNP): "txtListSNP"
% Here give the data file
% either: CNPSummary.mat, Output_deldup.mat or txt document
Indata = "20220623_LibBvalwS2d12_CNPSummary_wIndel";

% If you hand a .mat, where are your replicates of interest?
% If you want all samples, type []
dataset =  []; % [1:5]

% the plots and hotcolds are saved here
%savepath = "/home/isabel/sciebo/ResultsShared/DFE_HighThroughput/4_DNASeq/Promoters_Operons/";
savePath = "/media/isabel/IsabelFestplatte/3_DFEProject/4_DNASeq/Promoters_Operons";
ExpName = "LibSCBval";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Additional Files  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

donor = "Bval"; % "W23", "Bval", "Batro"

if strcmp(donor, "Bval")
    % Bval donor
    masterlist = '/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/mlBval2Bs166NCe_v1.txt';
    Accmm      = '/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/Bval_AccMM2Genes.mat';
elseif strcmp(donor, "W23")
    % W23 donor
    masterlist = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/mlW232Bs166NCe_v1.txt";
    Accmm = '/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/W23_AccMM2Genes.mat';
elseif strcmp(donor, "Batro")
    % BA donor
    masterlist = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/mlBatro2Bs166NCe_v1.txt";
    Accmm = '/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/Batro_AccMM2Genes.mat';
else
    error("Something went wrong with your donor declaration.. Is your donor really %s?", donor);
end


%%%%%%%     Load variables   %%%%%%%%%%%%
% which accmm genes to exclude? exclude_thr sets to which frac a gene must
% be accmm to be excluded; exclude_thr sets the lower limit
% % % if exclude_thr = 1, then only genes, that are complete accmm are exc.
% % % if exclude_thr = 1.1, then no gene is excluded
exclude_thr = 0.9;

% Which minimal number of MultiHits interests you? Set cutoff:
cutoff = 1;
recipsize = 4215607;

% % give me some information:

saveHotColdData = "OFF";
% are you interested in PROMOTER analysis??
PROM = "ON";


%% Here the data is loaded

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read in bedFile/ ml anc acc genome %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%Load recipient/donor specific master list
fid = fopen(masterlist); imp = textscan(fid,'%f %s %s');
fclose(fid);
refchr.pos=imp{1}; refchr.atcg=imp{3};
clear imp
clear prompt str fid ans imp

% % load the acc and mm genes for MultiHitStat
excAccmm = load(Accmm);
excAccmm = excAccmm.AccMM2Genes;

% load pseudoGenes

pseudoGenes = load(pseudoGenes);
pseudoGenes = pseudoGenes.pseudoGene;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read in main data of interest %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(IndataType,["Cdist","Adist","denovo","indel"])
    CNPSummary=load(IndataPath + Indata);
    if isempty(dataset)
        CNPSummary=CNPSummary.CNPSummary;
    else
        CNPSummary=CNPSummary.CNPSummary(dataset);
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
                Cluster{i,1}(1,:) = denovo{i,1}.pos(1,:)'; Cluster{i,1}(2,:) = denovo{i,1}.pos(1,:)';
                Cluster{i,2}(1,:) = denovo{i,1}.ref'; Cluster{i,3}(1,:) = denovo{i,1}.var';
               
            else
                Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
            end
        end

    elseif IndataType == "indel"
        for i=1:size(indel,1)
            if ~isempty(indel{i,1})
                Cluster{i,1}(1,:) = indel{i,1}.pos(1,:)'; Cluster{i,1}(2,:) = indel{i,1}.pos(1,:)';
                Cluster{i,2}(1,:) = indel{i,1}.ref'; Cluster{i,3}(1,:) = indel{i,1}.var';
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
    if isempty(dataset)
        DD.del=D.deldup.del;
        DD.dup=D.deldup.dup;
    else
        DD.del=D.deldup.del(dataset);
        DD.dup=D.deldup.dup(dataset);
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

elseif IndataType == "indel"
    load(IndataPath + Indata); %indelSummary loaded
    if ~isempty(dataset)
        indelSummary=indelSummary(dataset);
    end
    for i = 1 : length([indelSummary.pos])
        Cluster{i,1}(1,:) = indelSummary(i).pos; Cluster{i,1}(2,:) = indelSummary(i).pos;
        rep{i} = indelSummary(i).sample;
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
if ismember(IndataType,["del","dup","denovo"])
    exclude_thr = 1.1;
end

sampleNum = size(Cluster,1);

% for initializing structures
temp_string = strings(10000,1);
ts = num2cell(temp_string);
temp_noNum = (-1)* ones(10000,1);
tn = num2cell(temp_noNum);

%% HitOutput is generated

%Here for every replicate i all cluster c are checked for hits with genes and all different hits are written as entry into HitOutput.
%   -> This means: genes can appear multiple times per replicate -> this is
%      fixed in GeneOutput

HitOutput = struct('RepNo',tn,'Cluster',tn,'Type',ts,'Frac',tn,'locustag',ts,'Name',ts,'product',ts,'GeneDirection',ts,'Sample',ts);
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
            HitOutput(m).Type     = allTypes(t);
            HitOutput(m).RepNo    = i;
            HitOutput(m).Sample   = rep{i};
            HitOutput(m).Cluster  = c;
            HitOutput(m).Name     = string(bedFile(idxH).Name);
            HitOutput(m).locustag= string(bedFile(idxH).locus_tag);
            HitOutput(m).product  = string(bedFile(idxH).product);
            HitOutput(m).GeneDirection = string(bedFile(idxH).direction);
            
            if allTypes(t) == "tail"
                HitOutput(m).Frac = (bedFile(idxH).ende - Cluster{i}(1,c) + 1)/bedFile(idxH).L;
            end
            
            if allTypes(t) == "complete"
                HitOutput(m).Frac = 1;
            end
            
            if allTypes(t) == "in"
                HitOutput(m).Frac = (Cluster{i}(2,c) - Cluster{i}(1,c) + 1)/bedFile(idxH).L;
            end
            
            if allTypes(t) == "front"
                HitOutput(m).Frac = (Cluster{i}(2,c) - bedFile(idxH).start + 1)/bedFile(idxH).L;
            end
        end
        
        clear complete pin pfront ptail c co d pf pi pt
    end
end

% get rid of additional lines
HitOutput = HitOutput([HitOutput.locustag]~="");

%% GeneOutput and RepSummary are generated
% check for each hit gene if it also has a promoter that was hit!??

GeneOutput = struct('RepNo',tn,'Cluster',tn,'Type',ts,'Frac',tn,'locustag',ts,'Name',ts,'product',ts,...
    'GeneDirection',ts,'Sample',ts);
RepSummary = struct('RepNo',tn,'locustag_Hit',ts,'NameHit',ts,'Sample',ts);

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

% get rid of additional lines
GeneOutput = GeneOutput([GeneOutput.locustag]~="");
RepSummary = RepSummary([RepSummary.Sample]~="");


%% Analysis of PROMOTER!!

if PROM == "ON"
    
    % load the promoter regions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    operons = load(operonPath);
    operons = operons.operons_red;
    promoterLength = abs(operons(1).promoterEnde - operons(1).promoterStart) + 1;
    
    % in the following it will be interesting to know, which genes have
    % promoter informatoin at all ...
    % here: define list of all BSUs that have operon/promoter information
    genes_wPromInfo = unique([operons.BSUChain]);
    
    % % % first, HitOutput is generated for promoters
    %Here for every replicate i all cluster c are checked for hits with promoter and all different hits are written as entry into HitOutputProm.
    
    HitOutputProm = struct('RepNo',tn,'Cluster',tn,'Type',ts,'Frac',tn,'locustag',ts,'Sample',ts,'BSUChain',ts,'GeneChain',ts,'OperonDir',ts);
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
            if IndataType == "denovo"
                pin      = find((Cluster{i}(1,c) > [operons.promoterStart]) & (Cluster{i}(2,c) < [operons.promoterEnde]));
                
                allIdx   = [pin];
                allTypes = [repmat("in",1,numel(pin))];
            else
                ptail    = find(((Cluster{i}(1,c) > [operons.promoterStart]) & (Cluster{i}(1,c) < [operons.promoterEnde]) & (Cluster{i}(2,c) >= [operons.promoterEnde])));
                complete = find((Cluster{i}(1,c) <= [operons.promoterStart]) & (Cluster{i}(2,c) >= [operons.promoterEnde]));
                pin      = find((Cluster{i}(1,c) > [operons.promoterStart]) & (Cluster{i}(2,c) < [operons.promoterEnde]));
                pfront   = find((Cluster{i}(1,c) <= [operons.promoterStart]) & (Cluster{i}(2,c) > [operons.promoterStart]) & (Cluster{i}(2,c) < [operons.promoterEnde]));
                
                allIdx   = [ptail complete pin pfront];
                allTypes = [repmat("tail",1,numel(ptail)) repmat("complete",1,numel(complete))...
                    repmat("in",1,numel(pin)) repmat("front",1,numel(pfront))];
            end
            
            for t=1:numel(allIdx)
                idxH = allIdx(t);
                m = m+1;
                HitOutputProm(m).Type     = allTypes(t);
                HitOutputProm(m).RepNo    = i;
                HitOutputProm(m).Sample   = rep{i};
                HitOutputProm(m).Cluster  = c;
                HitOutputProm(m).BSUChain = string(operons(idxH).BSUChain);
                HitOutputProm(m).GeneChain= string(operons(idxH).GeneChain);
                HitOutputProm(m).locustag = string(operons(idxH).startGene);
                HitOutputProm(m).OperonDir= operons(idxH).direction;
                
                if allTypes(t) == "tail"
                    HitOutputProm(m).Frac = (operons(idxH).promoterEnde - Cluster{i}(1,c) + 1)/promoterLength;
                end
                
                if allTypes(t) == "complete"
                    HitOutputProm(m).Frac = 1;
                end
                
                if allTypes(t) == "in"
                    HitOutputProm(m).Frac = (Cluster{i}(2,c) - Cluster{i}(1,c) + 1)/promoterLength;
                end
                
                if allTypes(t) == "front"
                    HitOutputProm(m).Frac = (Cluster{i}(2,c) - operons(idxH).promoterStart + 1)/promoterLength;
                end
            end
            
            clear complete pin pfront ptail c co d pf pi pt
        end
    end
    
    
    % get rid of additional lines from initializing
    HitOutputProm = HitOutputProm([HitOutputProm.locustag]~="");
    
    % get a version of GeneOutput where the genes are connected with promoters
    GeneOutput_wProm = GeneOutput;
    
    for i=1:numel(GeneOutput_wProm)
        GeneOutput_wProm(i).Operon_BSUChain = {""};
        GeneOutput_wProm(i).OpInfo          = NaN;
        GeneOutput_wProm(i).Prom_Hit        = 0;
        GeneOutput_wProm(i).Prom_BSUChain   = "";
        GeneOutput_wProm(i).Prom_Dir        = "";
        GeneOutput_wProm(i).Prom_Cluster    = "";
        GeneOutput_wProm(i).Prom_Type       = "";
        GeneOutput_wProm(i).Prom_Frac       = "";
        
        % find the operons that the gene belons to
        operons_in_cells = {operons.BSUChain}; 
        idx_operon = find(cellfun(@(x) any(ismember(x,GeneOutput_wProm(i).locustag)),operons_in_cells));
        
        if ~isempty(idx_operon)
            GeneOutput_wProm(i).Operon_BSUChain = operons_in_cells(idx_operon);
            GeneOutput_wProm(i).OpInfo          = 1;
        else
            GeneOutput_wProm(i).Prom_Hit        = NaN;
        end
        
        % find the right replicate and the promoter that fits the gene
        idx_red = [HitOutputProm.RepNo] == GeneOutput_wProm(i).RepNo;
        BSUHit = cellfun(@(x) any(ismember(x,GeneOutput_wProm(i).locustag)),{HitOutputProm.BSUChain});
        
        % .. and put the information together!
        % promHit is in respect to HitOutputProm
        promHit = find(idx_red & BSUHit,true);
        
        if ~isempty(promHit)
            GeneOutput_wProm(i).Prom_Hit      = promHit;
            GeneOutput_wProm(i).Prom_Cluster  = HitOutputProm(promHit).Cluster;
            GeneOutput_wProm(i).Prom_Type     = HitOutputProm(promHit).Type;
            GeneOutput_wProm(i).Prom_BSUChain = HitOutputProm(promHit).BSUChain;
            GeneOutput_wProm(i).Prom_Dir      = HitOutputProm(promHit).OperonDir;
            GeneOutput_wProm(i).Prom_Frac     = HitOutputProm(promHit).Frac;
        end
        
        clear idx_red BSUHit promHit
    end
    
end

% how many hit genes have no operon/promoter info?
noInfo = sum(isnan([GeneOutput_wProm.OpInfo]));

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
        MultiHitStat(i).Samples   = GeneOutput(idx).Sample;
        
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
legend([h1 h2],{ExpName,ExpName + " wo AccGenes"})
