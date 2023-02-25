% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                    Covertgff3_to_bed.m                              % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% this script takes a gff/gff3 format file and extracts the fields of interest,
% which can be specified

%the way the data looks - EXAMPLE from .ggf3
% CP003909.1	Genbank	gene	101	1657	.	+	.	ID=gene-NGFG_00137;Name=NGFG_00137;gbkey=Gene;gene_biotype=protein_coding;locus_tag=NGFG_00137
% CP003909.1	Genbank	CDS	101	1657	.	+	0	ID=cds-EEZ46970.1;Parent=gene-NGFG_00137;Dbxref=NCBI_GP:EEZ46970.1;Name=EEZ46970.1;gbkey=CDS;locus_tag=NGFG_00137;product=chromosomal replication initiator protein dnaA;protein_id=EEZ46970.1;transl_table=11
% ...
% CP003909.1	Genbank	gene	13927	14012	.	-	.	ID=gene-NGFG_06000;Name=NGFG_06000;gbkey=Gene;gene_biotype=tRNA;locus_tag=NGFG_06000
% CP003909.1	Genbank	tRNA	13927	14012	.	-	.	ID=rna-NGFG_06000;Parent=gene-NGFG_06000;gbkey=tRNA;locus_tag=NGFG_06000;product=tRNA-Leu
% CP003909.1	Genbank	exon	13927	14012	.	-	.	ID=exon-NGFG_06000-1;Parent=rna-NGFG_06000;gbkey=tRNA;locus_tag=NGFG_06000;product=tRNA-Leu


clear all; close all
%%%%%%%%%%%%%%%%%%%%%%
% % specify the input file and the fields of interest
inPath = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/Ecoli/";
inDictName = "Ecoli";
inSuffix = ".gff3";
Header = 3;

outDictName = "Ecoli";

% if fasta was corrected with inserts/deletions, then shift all positions:
% if there is no shift, then
shift = zeros(5000000,1);

% ShiftPos = '/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/0_dictionaries/Bs166/CorrectFasta/ShiftBedPos_List.mat';
% shift = load(ShiftPos);
% shift = shift.shift;


% here, specify the fields that you are interested in
uniqList = struct('seqName',[],'start',[],'ende',[],'direction',[],'type1',[],'type2',[],'Name',[],'locus_tag',[],'product',[]);

% and where to put the output
outPath = inPath;
saveOutput = "ON";

fixWithSubtiWiki = "OFF";

%% Run the script
fid = fopen(inPath + inDictName + inSuffix); bed = textscan(fid,'%s %s %s %f %f %s %s %s %s','delimiter','\t','headerLines',Header);
fclose(fid);
genes.seqName=bed{1}; genes.type=bed{3}; genes.start=bed{4}; genes.ende=bed{5}; genes.direc=bed{7}; genes.info=bed{9};

searchFields = fieldnames(uniqList);
m=0;
for i=1:numel(genes.type)
    % the index m is changed only when an entry appears that belongs to a
    % new gene - most genes have 2-3 lines of information
    if i==1 || genes.start(i)~=genes.start(i-1)
        m = m+1;
        uniqList(m).seqName = string(genes.seqName{i});
        uniqList(m).start = genes.start(i) + shift(genes.start(i));
        uniqList(m).ende = genes.ende(i) + shift(genes.ende(i));
        uniqList(m).direction = string(genes.direc(i));
    end
    
    if strcmp(genes.type{i},'CDS')
        uniqList(m).type1 = string(genes.type(i));
    elseif strcmp(genes.type{i},'exon')
        uniqList(m).type1 = string(genes.type(i));
    elseif ~strcmp(genes.type{i},'exon') && ~strcmp(genes.type{i},'CDS')
        uniqList(m).type2 = string(genes.type(i));
    end
    
    if isempty(uniqList(m).type1)
        uniqList(m).type1 = "NA";
    end
    
    if isempty(uniqList(m).type2)
        uniqList(m).type2 = "NA";
    end
    
    % get the information out
    
    
    infocell = split(split(genes.info(i),";"),"=");
    for j=5:numel(searchFields)
        % first initialize
        
        % make an exception for the field GeneName, because there is no
        % name for each gene, and if not, take Name field from gene line
        if strcmp(searchFields(j),'Name')
            idx = find(strcmp(infocell,searchFields(j)));
            if ~isempty(idx) && (strcmp(genes.type(i),'gene') | strcmp(genes.type(i),'pseudogene'))

                    uniqList(m).(searchFields{j}) = string(infocell{idx,2});
              
            end
            
        else
            idx = find(strcmp(infocell,searchFields(j)));
            if ~isempty(idx)

                    uniqList(m).(searchFields{j}) = string(infocell{idx,2});
            
            end
        end
    end
    
    if isempty(uniqList(m).Name)
        uniqList(m).Name = "";
    end
    if isempty(uniqList(m).locus_tag)
        uniqList(m).locus_tag = "";
    end
end

% the problem can occur with subtilis, that the BSUs are not unqiue
% fix other problems
for i=1:numel([uniqList.start])
    if isempty(uniqList(i).locus_tag)
        uniqList(i).locus_tag = "";
        
    else
        %         spl = split(uniqList(i).old_locus_tag,"%");
        %
        %         if isempty(regexp(spl{1},"_"))
        %         spl2 = split(spl{1},"NGO");
        %         uniqList(i).locus_tag = "NGO" + "_" + string(spl2{2});
        %         else
        %          uniqList(i).locus_tag = string(spl{1});
        %         end
        uniqList(i).locus_tag = uniqList(i).locus_tag;
    end
end


%uniqList = rmfield(uniqList,"locus_tag")
%% test uniqueness of BSU, genename and positions

nonUniqBSU = cellfun(@(x) numel(find(strcmp(x,[uniqList.locus_tag]))),{uniqList(:).locus_tag});
num_nonUniqBSU = sum(nonUniqBSU > 1);
nonUniqGN  = cellfun(@(x) numel(find(strcmp(x,[uniqList.Name]))),{uniqList.Name});
num_nonUniqGN = sum(nonUniqGN > 1);

isListsorted = isequal(sort([uniqList(:).start]),[uniqList(:).start]);

fprintf('There are %i non-uniq BSU names and %i non-uniq genenames. \n', num_nonUniqBSU, num_nonUniqGN)

if isListsorted
    fprintf('Also, the uniqList is sorted according to the start position. \n')
else
    fprintf('Also, the uniqList is NOT sorted according to the start position. \n')
end

maskNoLocusTag = [uniqList.locus_tag] == "";

fprintf('There are %i entries without a locustag that are being deleted. \n', sum(maskNoLocusTag));

uniqList = uniqList(~maskNoLocusTag);

%% for subtilis, read in the subti wiki 

if fixWithSubtiWiki == "ON"
% fix the problem that parts of the annotation are wrong
% find the gene names that appear more often than once and check their BSU
% - genename equivalence!

fid = fopen("/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/0_dictionaries/Bs166_unchanged/" + "Bs166Names_from_SubtiWiki_June2021.csv"); 
subtiWikiInfo = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','delimiter',',','headerLines',1);
fclose(fid);
subti.BSU=subtiWikiInfo{2}; subti.GN=subtiWikiInfo{3}; 

idx_fix = find(nonUniqGN>1);

uniqList_fix = uniqList;

for i=1:numel(idx_fix)
   % this is the problematic entry in the uniqList 
   idx  = idx_fix(i);
   lt   = uniqList_fix(idx).locus_tag;
   gene = uniqList_fix(idx).Name; 
    
   % find where this BSU is in the subti list
   test = cellfun(@(x) find(strcmp(lt,string(x))),{subti.BSU},'UniformOutput',false);
   
   if ~isempty(test{1})
   idx_inSubti = cellfun(@(x) find(strcmp(lt,string(x))),{subti.BSU}); 
   gene_inSubti = string(subti.GN{idx_inSubti});
   
   uniqList_fix(idx).Name = gene_inSubti;
   
   fprintf('Gene %s, %s is now %s, according to subtiWiki. \n', gene, lt, gene_inSubti)
   
   elseif isempty(test{1})
       
       fprintf('No equivalence for %s according to subtiWiki. \n', gene)
   end
    
end
uniqList = uniqList_fix;
end



%% write to file

if saveOutput == "ON"
    % Option: filter for specific type, exon or CDS
    maskCDS   = [uniqList.type1]=="CDS";
    maskExon  = [uniqList.type1]=="exon";
    maskOther = [uniqList.type1]~="CDS" & [uniqList.type1]~="exon";
    
    uniqList_CDS = uniqList(maskCDS);
    uniqList_exon = uniqList(maskExon);
    uniqList_other = uniqList(maskOther);
    
    % save as mat -- the whole list
    save(outPath + outDictName + ".bed.mat",'uniqList')
    % as table in .txt
    T=struct2table(uniqList);
    writetable(T,outPath + outDictName + ".bed.txt",'Delimiter','\t')
    
    % save as mat -- only CDS
    save(outPath + outDictName + "_onlyCDS.bed.mat",'uniqList_CDS')
    % as table in .txt
    T_CDS=struct2table(uniqList_CDS);
    writetable(T_CDS,outPath + outDictName + "_onlyCDS.bed.txt",'Delimiter','\t')
    
    % save as mat -- only exon
    save(outPath + outDictName + "_onlyExon.bed.mat",'uniqList_exon')
    % as table in .txt
    T_exon=struct2table(uniqList_exon);
    writetable(T_exon,outPath + outDictName + "_onlyExon.bed.txt",'Delimiter','\t')
    
end
