% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                    Annotate_with_otherSpecies.m                              % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% With this script you can take one badly annotated species (originally created fpr Ngo MS11)
% and annotate it with different other species.

% Input:
% -- a list of blast results, where the badly annotated species is used as a Query and blasted against
%    the other species (Subject) gene by gene
% -- bed file of each species that should be used for the annotation


%the way the data looks - EXAMPLE
% Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
%NGFG_00137	NC_003112.2	95.568	1557	69	0	1	1557	2005945	2004389	0.0	2494
%NGFG_00138	NC_003112.2	96.739	1104	36	0	1	1104	2004154	2003051	0.0	1840
%NGFG_00139	NC_003112.2	94.275	2061	112	4	1	2058	2002040	1999983	0.0	3147
%NGFG_00140	NC_003112.2	92.568	444	33	0	1	444	1999814	1999371	0.0	638
%pseudogene5665	NC_003112.2	97.099	517	14	1	1	517	1999364	1998849	0.0	870
%NGFG_00143	NC_003112.2	95.606	2640	110	3	1	2637	1998585	1995949	0.0	4228
%pseudogene9307	NC_003112.2	96.469	793	26	2	1	793	1995815	1995025	0.0	1308
%...
clear all; close all

%% specify the input file and the fields of interest
% in blast, the query Ref is blasted against the subject Sib (sibling1)
% and possible also other species -> this gives a blast list and for both
% the Ref and Sib there are bed files
basePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/";
blastList = basePath + "1B_MS11_ncbi/Ann_Ngo2Sibs/singleGenesBlast2Sibs/" + ["Blast_NgonGenes_to_Nmen_editted.tab",...
    "Blast_NgonGenes_to_FA1090_editted.tab",...
    "Blast_MS11NCBI_singleGenes_nucl_2_NZ_AP0230691_editted.tab"];

BedFileQuery = basePath + "/1A_MS11_ManualLists/MS11.bed.mat";

BedFileSubject = basePath + ["4_Nmeningitidis_MC58/ManualLists/Nmeningitidis.bed.mat",...
    "3_FA1090/ManualLists/FA1090.bed.mat",...
    "5_NZ_AP0230691_MS11_2022/ManualLists/NZ_AP0230691.bed"];

Query = "Ngo";

Subject = ["Nmen" "FA1090" "NgoMS11_new"];
chrQuery = [2233640];
chrSubject = [2272360,2153922,2167602];

MultiMapMode = "OFF"; % this you have to set, because the function B_SortBlstHits takes this variable!
FixMode      = "OFF";

% these types are excluded according to the query and subject
excType_inQuery = "NA";
excType_inSubject = "NA";

% hits that are under a length fraction of this for Query and Subject are
% excluded
excLfracQuery   = 0.5;
excLfracSubject = 0.5;

% and where to put the output
outPath = basePath + "1A_MS11_ManualLists/";
SaveOutput = "OFF";

%% Run the script SortBlastHits.m
blastResult = struct("Query2",[],"out",[]);
for i=1:numel(blastList)
    blastResult(i).Query2 = Subject(i);
    [blastMatch,blastMatch_reduced] = B_SortBlastHits(blastList(i),BedFileQuery,BedFileSubject(i),chrSubject(i),MultiMapMode,FixMode,excType_inQuery,excType_inSubject,excLfracQuery,excLfracSubject);
    blastResult(i).out = blastMatch_reduced;

end

%% Nooooowwww the problem is that we want to get a 
% list in which all genes of the first species are listed!!
% not just the ones that had blast hits ..
% So now, we go through the bedQ list and write collectList, in which we
% add the information from blastMatch_reduced, if available

bedQuery = load(BedFileQuery);
bedQuery = bedQuery.uniqList;

temp_string = strings(numel([bedQuery.seqName]),1);
ts = num2cell(temp_string);

temp_NAN    = nan(numel([bedQuery.seqName]),1);
tn = num2cell(temp_NAN);

% initialize the collectList structure, so that a variable number of other
% specie's annotations can be used
Fields = ["locustag_in_","GN_in_","type_in_","product_in_"];

m = 0;
initStruct = struct('collect',[]);

%fist initialize the query fields
for i=1:4
    m = m +1;
    initStruct(m).collect = Fields(i) + Query;
    m = m +1;
    initStruct(end+1).collect = ts;
end
initStruct(end+1).collect = "BlastHit";
initStruct(end+1).collect = [];

% then initialize the other strains
for s = 1:numel(chrSubject)
    for i=1:4
        initStruct(end+1).collect = Fields(i)   + Subject(s);
        initStruct(end+1).collect = ts;
    end
end

collectList = struct(initStruct(:).collect);
bedFields = ["locus_tag","Name","type2","product"];
QuFields = Fields + Query;

for i=1:numel(bedQuery)
    
    % first get the information from the query bed file
    for j=1:numel(QuFields)
        if ~isempty(bedQuery(i).(bedFields(j)))
        collectList(i).(QuFields(j)) = bedQuery(i).(bedFields(j));
        end
    end    
    
    for s = 1:numel(Subject)
        
        blastMatch_temp = blastResult(s).out;
        
        idx_inBlast = find(strcmp(collectList(i).(QuFields(1)),[blastMatch_temp.locustag_in_Query]));
        
        if ~isempty(idx_inBlast)
            
            collectList(i).BlastHit = [collectList(i).BlastHit Subject(s)];
            
            collectList(i).("type_in_"     + Subject(s)) = blastMatch_temp(idx_inBlast).type_in_Subject;
            collectList(i).("product_in_"  + Subject(s)) = blastMatch_temp(idx_inBlast).product_in_Subject;
            collectList(i).("locustag_in_" + Subject(s)) = blastMatch_temp(idx_inBlast).locustag_in_Subject;
            collectList(i).("GN_in_" + Subject(s)) = blastMatch_temp(idx_inBlast).GN_in_Subject;
            
        end
    end
          
end

%% Now as a last step give, for query and each subject, the columns GN_in_... and product_in_... to  
% function called: Find_GeneNames.m, that will extract, if possible, the gene names
goThrough = [Query Subject];

for i=1:numel(goThrough)
GN  = [collectList.("GN_in_"      + goThrough(i))];
Pro   = [collectList.("product_in_" + goThrough(i))];
excWords = ["beta","II","RNA","rrf","M23","L" + [1:40],"S" + [1:40],"IF-1","IF-3","IF-2","L38","L28","ill","iii","type","SMC","HII",...
    "[Fe]","B561","2/3","5/6)","p.1B","h.8","h-NS"];

out = N2_Find_GeneNames(GN,Pro,excWords);

% build it into the collectList 
[collectList.("GN_in_"+ goThrough(i))] = deal(out.gene);
end

%% additional stuff for Ngo (14.10.2022)
% you can delete this!!!
collectList_onlyGenes = rmfield(collectList,{'product_in_FA1090','product_in_Ms11_2022','product_in_Ngo','product_in_Nmen','BlastHit','locustag_in_FA1090','locustag_in_Ms11_2022','locustag_in_Nmen'});

%%
% % save as mat
if SaveOutput == "ON"

    save([outPath + "NgoBlast2_" + Subject(1) + "AND" + Subject(2) + ".mat"],'collectList')
    % as table in .txt
    T=struct2table(collectList);
    writetable(T,[outPath + "NgoBlast2_" + Subject(1) + "AND" + Subject(2) + ".csv"],'Delimiter','\t')

    % only wo CDS
    maskCDS =  [collectList.type_in_Ngo]=="gene";
    collectList_CDS = collectList(maskCDS);

    collectList_CDS = rmfield(collectList_CDS,{'BlastHit','type_in_Nmen','product_in_Nmen',...
        'type_in_FA1090','product_in_FA1090'});

    save([outPath + "NgoBlast2_" + Subject(1) + "AND" + Subject(2) + "_onlyCDS.mat"],'collectList_CDS')
    % as table in .txt
    T_CDS=struct2table(collectList_CDS);
    writetable(T_CDS,[outPath + "NgoBlast2_" + Subject(1) + "AND" + Subject(2) + "_onlyCDS.csv"],'Delimiter','\t')

end