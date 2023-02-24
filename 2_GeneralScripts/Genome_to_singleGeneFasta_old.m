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

clear all; close all

GenomeFasta = '/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/BA1942.fasta';
BedFile     = '/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/Batro_June2021.bed.mat';
outPath     = "/home/isabel/Desktop/";

species = "Bvallismortis";

% load data
indata = fastaread(GenomeFasta);
nt = indata.Sequence;

bed = load(BedFile);
bed = bed.uniqList;

printFigs = "ON";

%% here the nt sequences are written to the single gene output
% each gene is written as separate contig with line width of 70

if exist([outPath + species + "_singleGenes_nucl.fasta"], 'file')==2
  delete([outPath + species + "_singleGenes_nucl.fasta"]);
end

fileID = fopen([outPath + species + "_singleGenes_nucl.fasta"],'a');

for i=1:numel([bed(:).start])
    
    % write which gene you are looking at
    if isempty(bed(i).Name) || isempty(bed(i).product)
    fprintf(fileID,'%s \n',">" + string(bed(i).locus_tag));      
    else
    fprintf(fileID,'%s \n',">" + string(bed(i).locus_tag) + " " + string(bed(i).Name) + " " + string(bed(i).product));      
    end
    % get the nt sequence:
        % check in which direction the gene is read out
        % if (-), then make a reverse complimentary
    if strcmp(bed(i).direction,'+')
        sequence = nt(bed(i).start:bed(i).ende);
    elseif strcmp(bed(i).direction,'-')
        sequence = seqrcomplement(nt(bed(i).start:bed(i).ende));
    end
        
    width = 70; % width of the lines in the file
    check  = [];
    check2 = [];
    for j=1:ceil(numel(sequence)/width)
              
        if j*width <= numel(sequence)
            seq = sequence(1+width*(j-1):j*width);
            
            % in check and check2 you can see what you wrote to file
            check(j,:) = 1+width*(j-1):j*width;
            check2(j,:)= sequence(1+width*(j-1):j*width);
            
        elseif  j*width > numel(sequence)
            seq = sequence(1+width*(j-1):end);
        end
        
        fprintf(fileID,[repmat('%c',1,numel(seq)),'\n'],seq);
    end
end

fclose(fileID);

%% here the nt sequences are translated to amino acid sequences and 
% written to output

if exist([outPath + species + "_singleGenes_aa.fasta"], 'file')==2
  delete([outPath + species + "_singleGenes_aa.fasta"]);
end

fileID = fopen([outPath + species + "_singleGenes_aa.fasta"],'a');

% collect the information on the protein integrity: does sequence start
% with M and end with*? what is its AA sequence?
AAusage = struct('Name',num2cell(strings(numel(bed),1)),'locus_tag',num2cell(strings(numel(bed),1)),'type1',num2cell(strings(numel(bed),1)),...
    'type2',num2cell(strings(numel(bed),1)),'StartAA',[],'StopAA',[],'SequenceAA',[]);
AAfields = fieldnames(AAusage);

for i=1:numel([bed(:).start])
    
    % write which gene you are looking at
    fprintf(fileID,'%s\n',">" + string(bed(i).locus_tag) + " " + bed(i).Name + " " + string(bed(i).product));  
    
    for k=1:4
        if ~isempty(bed(i).(AAfields{k}))
    AAusage(i).(AAfields{k}) = string(bed(i).(AAfields{k}));
        end
    end
    
    % now check in which direction the gene is read out
    if strcmp(bed(i).direction,'+')
        sequence_tmp = nt(bed(i).start:bed(i).ende);
        sequence = nt2aa(sequence_tmp, 'ACGTOnly','false');
    elseif strcmp(bed(i).direction,'-')
        sequence_tmp = seqrcomplement(nt(bed(i).start:bed(i).ende));
        sequence = nt2aa(sequence_tmp,'ACGTOnly','false');
    end

        AAusage(i).StartAA    = string(sequence(1)); 
        AAusage(i).StopAA     = string(sequence(end));
        AAusage(i).SequenceAA = sequence;
    
    width = 70;
    check  = [];
    check2 = [];
    for j=1:ceil(numel(sequence)/width)
              
        if j*width <= numel(sequence)
            seq = sequence(1+width*(j-1):j*width);
            check(j,:)=1+width*(j-1):j*width;
            check2(j,:)=sequence(1+width*(j-1):j*width);
            
        elseif  j*width > numel(sequence)
            seq = sequence(1+width*(j-1):end);

        end
        
        fprintf(fileID,[repmat('%c',1,numel(seq)),'\n'],seq);
    end 
    
end

fclose(fileID);

%% analyse all amino acids and their usage/ start and stop amino acids

% create masks for the type of gene that you are interested in
onlyCDS  = [AAusage.type1] == "CDS";
onlygene = [AAusage.type2] == "gene";
RNA      = contains([AAusage.type2],'RNA');
noEntry  = [AAusage.type1] == "";

%for which mask do you want to look at the usage of codons?
mask = onlyCDS & onlygene;

% plot the usage of all amino acids
figure(1)
title(['Usage of all amino acids in CDS genes of ' [' '] species], 'Interpreter', 'none'); hold on
[a,b,c] = unique(strcat(AAusage(mask).SequenceAA));
for i=1:length(a)
    uniqAA{i} = a(i);
end
figure(1)
h=histogram(c);
xticks(1:length(uniqAA))
xticklabels(uniqAA);
set(gca,'Yscale','lin')
ylim([0 200000])

% plot the usage of start amino acids
figure(10)
title(['Usage of start amino acids in CDS genes of ' [' '] species], 'Interpreter', 'none'); hold on
[a,b,c] = unique([AAusage(mask).StartAA]);
for i=1:length(a)
    uniqStart{i} = a(i);
end
figure(10)
h=histogram(c);
xticks(1:length(uniqStart))
xticklabels(uniqStart);
set(gca,'Yscale','log')
clear uniqStart

figure(11)
title(['Usage of stop amino acids in CDS genes of ' [' '] species], 'Interpreter', 'none'); hold on
[d,f,g] = unique([AAusage(mask).StopAA]);
for i=1:length(d)
    uniqStop{i} = d(i);
end
figure(11)
h=histogram(g);
xticks(1:length(uniqStop))
xticklabels(uniqStop);
set(gca,'Yscale','log')
clear uniqStart

%%
% plot the figures

if printFigs == "ON"
    print(figure(1),'-painters','-dpng', outPath + "AllAA_Usage_" + species + "_all.png")
    print(figure(10),'-painters','-dpng', outPath + "StartAA_Usage_" + species + "_all.png")
    print(figure(11),'-painters','-dpng', outPath + "StopAA_Usage_" + species + "_all.png")
end
