%% the idea of this script is to ...
%
%
% in: 
% % %  -- list of operons from SubtiWiki
% % %  -- bed file of the species in question
% % %  -- a name with the Names connected to BSU and synonyms
% 
% out: 
% % %  -- list of operons with start end and promoter region

close all; clear all
%% set paths and read in data
basePath = "/home/isabel/sciebo/ResultsShared/kleinesPaper/allLists/";

promoterPath = basePath + "20210825_PromoterRegions.mat";
recipbed       = basePath + "allLists/" + "Bs166NCe_June2021.bed.mat";
% bedFile     = basePath + "Bs166NCe_June2021.bed.mat";
% RefFasta    = basePath + "Bs166NCe.fasta";

dataInPath = "/home/isabel/sciebo/ResultsShared/DFE_HighThroughput/4_DNASeq/Promoters_Operons/";

dataset(1) = "20210902_GeneOutput_wProm_LibSCBval.mat";
%dataset(2) = "20210901_GeneOutput_wProm_LibSCW23Outlier.mat";

 % fid = fopen(Names2BSU_path);
% Names2BSU = readtable(Names2BSU_path);
% fclose(fid); clear fid
% 
% bedFile = load(bedFile);
% bedFile = bedFile.uniqList;
% 
% fasta = fastaread(RefFasta);
% fasta = fasta.Sequence;

recipsize = 4215607;

%Samples = ["LibSCW23" + string([1:10]), "LibSCW2318","LibSCW2332","LibSCW2387"];
Samples = ["LibSCBval" + [1 4 5 7 9 17 89 92]];

ExpName = "LibSCBval";

savePath = "/home/isabel/sciebo/ResultsShared/DFE_HighThroughput/4_DNASeq/Promoters_Operons/plots/";

%% load data

% load the gene data
for i=1:numel(dataset)
    A = load(dataInPath + dataset(i)) ;
    if i==1
        GeneOutcollect = A.GeneOutput_wProm;
    else
        GeneOutcollect = [GeneOutcollect; A.GeneOutput_wProm];
    end
end


% load the promoter regions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
promoter = load(promoterPath);
promoter = promoter.operons_red;
promoterLength = abs(promoter(1).promoterEnde - promoter(1).promoterStart) + 1;

% in the following it will be interesting to know, which genes have
% promoter informatoin at all ...
genes_wPromInfo = unique([promoter.BSUChain]);


%% 

% how many genes have a promoter that is hit??

% replicate = [GeneOutput_wProm.RepNo];

genesHit = zeros(numel(Samples),1);
for i=1:numel(Samples)
    clear rep_mask geneHit_wPromInfo genesHit_wProm
    
    smp_mask          = [GeneOutcollect.Sample] == Samples(i);
    GenesHere         = GeneOutcollect(smp_mask);
    
    % all genes that were hit
    genesHit(i)       = numel(GenesHere);
    
    if genesHit(i) == 0
        geneHit_wPromInfo = 0; geneHit_woProm = 0;
    else
        % distinguish between those with info and wo info
        geneHit_wPromInfo_mask  = ismember([GenesHere.locustag],genes_wPromInfo);
        geneHit_wPromInfo       = sum(geneHit_wPromInfo_mask);
        geneHit_woPromInfo      = genesHit(i)-geneHit_wPromInfo;
        
        % here only look at those, that have information
        geneHit_wProm   = sum([GenesHere(geneHit_wPromInfo_mask).Prom_Hit]~=0);
        geneHit_woProm  = sum([GenesHere(geneHit_wPromInfo_mask).Prom_Hit]==0);
    end
    bar_collect(i,:)      = [geneHit_wProm,geneHit_woProm,geneHit_woPromInfo];
    bar_collect_norm(i,:) = 100 * [geneHit_wProm,geneHit_woProm,geneHit_woPromInfo]/genesHit(i);

end

cmp(1,:) = [0 68 136]/255;
cmp(2,:) = [221 170 51]/255;
cmp(3,:) = [187 85 102]/255;

figure(102); hold on; set(gcf, 'Position', [550 350 600 750], 'Renderer', 'painters');

subplot(2,1,1)

bar([1:numel(Samples)],genesHit,'FaceColor','none','EdgeColor',[0.4 0.4 0.4],'LineWidth',1.4)
set(gca,'XTick',[])
ylim([0.5 45])
ylabel('by replacement in total')
%legend({'gene and promoter hit','only gene hit'})

subplot(2,1,2)

b = bar(bar_collect_norm,'stacked');
b(1).FaceColor = cmp(1,:);
b(2).FaceColor = cmp(2,:);
b(3).FaceColor = cmp(3,:);
%set(gca,'XTick',[1:numel(Samples)],'XTickLabel',cellstr(Samples),'XTickLabelRotation',45);
set(gca,'XTick',[1:numel(Samples)],'XTickLabel',{'LibSCBval1','LibSCBval37','LibSCBval49','LibSCBval7','LibSCBval9','LibSCBval17','LibSCBval89','LibSCBval92'},'XTickLabelRotation',45);
ylim([0 102]);
legend({'including promoter','without','no info on promoter'},'Location','northoutside','NumColumns',3)

subplot(2,1,1)
title("Genes that are hit in " + ExpName + "..")