% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                   Genome_to_singleGeneFasta.m             % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%REQUIRED: bioinformatic toolbox

% This script takes as
%
% INPUT:
% ------ a whole genoma fasta and a bed file (in .mat) of one species (e.g. Bacillus subtilis)
%
% It then takes the start and end position of each gene from the bed file and retrieves it's
% nucleotide sequence from the fasta. These are written as
%
% OUTPUT:
%  ----- to a .fasta file of the form
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   >BSU_00010 dnaA chromosomal replication initiator informational ATPase
%   ATGGAAAATATATTAGACCTGTGGAACCAAGCCCTTGCTCAAATCGAAAAAAAGTTGAGCAAACCGAGTT
%   TTGAGACTTGGATGAAGTCAACCAAAGCCCACTCACTGCAAGGCGATACATTAACAATCACGGCTCCCAA
%   .....
%   >BSU_00020 dnaN DNA polymerase III (beta subunit)
%   ATGAAATTCACGATTCAAAAAGATCGTCTTGTTGAAAGTGTCCAAGATGTATTAAAAGCAGTTTCATCCA
%   ....
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  and
%  ---- to it converts the nucleotides to amino acids and writes it to a
%  .fasta
%   --> this can then be e.g. blasted against a genome
%
% Further stuff:
%
% keep in mind: if the program finds that a gene is written in the other
% direction, it automatically creates the reverse complementary
%
addpath("/home/isabel/Documents/Doktorarbeit_Mai2022/SCRIPTS/matlab2TikZ/src/")
clear all; close all

% GenomeFasta = 'Bs166NCe.fasta';
% BedFile     = '/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/';
% GenomeFasta = '/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/';
% BedFile     = ['/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/W23wt_June2021.bed.mat'];
basePathBac = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/";
basePathNgo = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/";
GenomeFasta = [basePathBac + ["Bs166NCe.fasta" "W23wt.fasta" "Bvallismortis.fasta"] basePathNgo + "1B_MS11_ncbi/MS11.fasta"];
BedFile     = [basePathBac + ["Bs166NCe_June2021.bed.mat" "W23wt_June2021.bed.mat" "Bvallismortis_June2021.bed.mat"] basePathNgo + "1A_MS11_ManualLists/MS11.bed.mat"];

%TrafoExp   = ["Wns","Vns","BAns","Geons"];
DonorColor = [0 0 0; 51 34 136; 136 34 85; 221 221 221]/265;
gray = [221 221 221]/265;

species = ["Bsub" "Bspiz" "Bval" "MS11"];
printFigs = "OFF";

%% now collect the data of how often which aa appears at the start and end!!
AAStart = ["M" "V" "I" "other"];
AAStop  = ["*" "other"];

for s=1:numel(species)

    % load data
    indata = fastaread(GenomeFasta(s));
    nt{s} = indata.Sequence;
    clear indata

    bed_temp = load(BedFile(s));
    collectBed(s).bed = bed_temp.uniqList;
    collectBed(s).species = species(s);
    clear bed_temp



    % collect the information on the protein integrity: does sequence start
    % with M and end with*? what is its AA sequence?

    % go through the samples!
    collect(s).species = species(s);
    collect(s).AAusage = struct('Name',num2cell(strings(numel(collectBed(s).bed),1)),'locus_tag',num2cell(strings(numel(collectBed(s).bed),1)),'type1',num2cell(strings(numel(collectBed(s).bed),1)),...
        'type2',num2cell(strings(numel(collectBed(s).bed),1)),'StartAA',[],'StopAA',[],'SequenceAA',[]);
    AAfields = fieldnames(collect(s).AAusage);

    for i=1:numel([collectBed(s).bed.start])

        for k=1:4
            if ~isempty(collectBed(s).bed(i).(AAfields{k}))
                collect(s).AAusage(i).(AAfields{k}) = string(collectBed(s).bed(i).(AAfields{k}));
            end
        end

        % now check in which direction the gene is read out
        if strcmp(collectBed(s).bed(i).direction,'+')
            sequence_tmp = nt{s}(collectBed(s).bed(i).start:collectBed(s).bed(i).ende);
            sequence = nt2aa(sequence_tmp, 'ACGTOnly','false');
        elseif strcmp(collectBed(s).bed(i).direction,'-')
            sequence_tmp = seqrcomplement(nt{s}(collectBed(s).bed(i).start:collectBed(s).bed(i).ende));
            sequence = nt2aa(sequence_tmp,'ACGTOnly','false');
        end

        collect(s).AAusage(i).StartAA    = string(sequence(1));
        collect(s).AAusage(i).StopAA     = string(sequence(end));
        collect(s).AAusage(i).SequenceAA = sequence;

    end

    % analyse all amino acids and their usage/ start and stop amino acids

    % create masks for the type of gene that you are interested in
    onlyCDS  = [collect(s).AAusage.type1] == "CDS";
    onlygene = [collect(s).AAusage.type2] == "gene";
    RNA      = contains([collect(s).AAusage.type2],'RNA');
    noEntry  = [collect(s).AAusage.type1] == "";

    %for which mask do you want to look at the usage of codons?
    mask{s}            = onlyCDS & onlygene;
    collect(s).NumCDS  = numel([collect(s).AAusage(mask{s})]);


    for i=1:numel(AAStart)-1
        collect(s).StartAA(i) = sum([collect(s).AAusage(mask{s}).StartAA]==AAStart(i));
    end
    collect(s).StartAA(4)      = numel([collect(s).AAusage(mask{s})]) - sum(collect(s).StartAA(1:3));
    collect(s).StartAA_perc = collect(s).StartAA/collect(s).NumCDS;

    for i=1:numel(AAStop)-1
        collect(s).StopAA(i) = sum([collect(s).AAusage(mask{s}).StopAA]==AAStop(i));
    end
    collect(s).StopAA(2)      = numel([collect(s).AAusage(mask{s})]) - sum(collect(s).StopAA(1:1));
    collect(s).StopAA_perc = collect(s).StopAA/collect(s).NumCDS;

end

%% plot the result!!!
StartAA_perc = vertcat(collect.StartAA_perc);
StopAA_perc  = vertcat(collect.StopAA_perc);

figure(1); hold on
b=bar(StartAA_perc')
set(gca,'XTick',[1:numel(AAStart)],'XTickLabel',AAStart)
legend(species,'Interpreter','none')
for i=1:size(StartAA_perc',2)
b(i).FaceColor =DonorColor(i,:);
end
ylabel('percentage of start codons')


figure(2); hold on
bb=bar(StopAA_perc')
set(gca,'XTick',[1:numel(AAStop)],'XTickLabel',AAStop)
legend(species,'Interpreter','none')
for i=1:size(StopAA_perc',2)
bb(i).FaceColor =DonorColor(i,:);
end
ylabel('percentage of stop codons')


