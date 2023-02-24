% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                    BlastHits_2_Lists.m                              % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%%%% This script takes a list of blast results in table form and the bed file of the two
%%%% species (query & subject) and converts it to a blastMatch list.
%%%%%% The script can be run for blast lists between different species (e.g. gene equivalence 
%%%%%% search between species) or between the same species (MULTIMAPPER)

%%%%%%%% What the script does %%%%%%%
%%%%  -- Part1: the blast hit list is sorted

%the way the data looks - EXAMPLE 
% Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
% BSU_00010	NZ_CP026362.1	95.007	1342	65	2	1	1341	3580217	3578877	0.0	2106
% BSU_00020	NZ_CP026362.1	91.909	1137	92	0	1	1137	3578690	3577554	0.0	1591
% BSU_00030	NZ_CP026362.1	94.907	216	11	0	1	216	3577422	3577207	7.17e-94	339
% BSU_00040	NZ_CP026362.1	94.519	1113	61	0	1	1113	3577191	3576079	0.0	1718
% BSU_00050	NZ_CP026362.1	97.967	246	5	0	1	246	3576061	3575816	1.72e-120	427
% BSU_00060	NZ_CP026362.1	93.166	1917	131	0	1	1917	3575761	3573845	0.0	2815
% BSU_00060	NZ_CP026362.1	80.556	108	19	2	280	386	1504401	1504295	1.30e-15	82.4
% BSU_00070	NZ_CP026362.1	92.660	2466	181	0	1	2466	3573634	3571169	0.0	3552
% BSU_rRNA_1	NZ_CP026362.1	99.743	1554	4	0	2	1555	3484059	3482506	0.0	2848
% BSU_rRNA_1	NZ_CP026362.1	99.678	1554	5	0	2	1555	3489986	3488433	0.0	2843
% BSU_rRNA_1	NZ_CP026362.1	99.678	1554	4	1	2	1555	2553671	2552119	0.0	2841
%...

clear all; close all

%% specify the input file and the fields of interest
% in blast, the query is blasted against the subject
% and possible also other species -> this gives a blast list and for both
% the Query and Subject there are bed files
basePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/1A_MS11_ManualLists/";
blastList      = basePath + "MULTIMAPPER_Ms11toMs11/Blast_MS11Genes_2_MS11_editted.tab";
BedFileQuery   =  basePath +"MS11.bed.mat";
BedFileSubject =  basePath +"MS11.bed.mat";

Query = 'MS11';
Subject = 'MS11';
 
MultiMapMode = "ON";

% these types are excluded according to the query and subject
excType_inQuery = ["pseudogene","RNA"];
excType_inSubject = 'NA';

% hits that are under a length fraction of this for Query and Subject are
% excluded
excLfracQuery   = 0.3;
excLfracSubject = 0.3;

chrSubject = 4286362;

saveOutput = "OFF";

% and where to put the output
outPath = '/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/MS11_ncbi_IR/MULTIMAPPER_Ms11toMs11/';

%% Run the SortBlastHits.m function 
FixMode = "OFF";
[blastMatch,blastMatch_reduced,blastMatch_mm] = B_SortBlastHits(blastList,BedFileQuery,BedFileSubject,chrSubject,MultiMapMode,FixMode,excType_inQuery,excType_inSubject,excLfracQuery,excLfracSubject);

%% 
if MultiMapMode == "ON"
    ident_MM = 0.99; % this identity defines a mm -- it is relatively strict, 
                     % because the Lfrac defining a mm is not so strict
    for i=1:numel(blastMatch_mm)
        
            identCheck = [blastMatch_mm(i).otherHits(:).ident] >= ident_MM;
        
            blastMatch_mm(i).otherHits          = blastMatch_mm(i).otherHits(identCheck);
            blastMatch_mm(i).otherHitsNum       = sum(identCheck);
          
       
    end
end

blastMatch_mm = blastMatch_mm([blastMatch_mm.otherHitsNum] > 0);
blastMatch_mm = blastMatch_mm(~contains([blastMatch_mm.type_in_Query],"RNA"));
%% save the outputs

if saveOutput == "ON"
    if MultiMapMode == "ON"
        save(outPath + "MultiMapper_" + Query + "_2_" + Subject + ".mat",'blastMatch_mm')
        T_MM=struct2table(blastMatch_mm);
        writetable(T_MM,[outPath 'MultiMapper_' Query '_2_' Subject '.csv'],'Delimiter','\t')
    else
        % save as mat
        save([outPath 'GeneEquivalence_' Query '_' Subject '.mat'],'blastMatch_reduced')
        % as table in .txt
        T=struct2table(blastMatch_reduced);
        writetable(T,[outPath 'GeneEquivalence_' Query '_' Subject '.csv'],'Delimiter','\t')
       
    end
    
end


    %% Plot the distributions
if MultiMapMode == "OFF"    
    % Distribtuion of identities
    figure(1)
    h=histogram([blastMatch_reduced(:).ident],'FaceColor','none'); hold on
    
    plot([mean([blastMatch_reduced(:).ident]) mean([blastMatch_reduced(:).ident])],[0 300],'LineWidth',3)
    
    % Look at the lengths of the genes
    LengthQuery = [blastMatch_reduced(:).endeQuery]-[blastMatch_reduced(:).startQuery]+1;
    LengthSubject = abs([blastMatch_reduced(:).endeSubject]-[blastMatch_reduced(:).startSubject])+1;
    
    figure(2)
    scatter(LengthQuery,LengthSubject)
    
    % look at how length and identity go hand in hand
    figure(3)
    scatter([blastMatch_reduced(:).ident],LengthQuery,'Marker','+','MarkerEdgeColor',[0.7 0.7 0.7])
    
    %% Look at different subsets
    numSam = 10000;
    for i=1:numSam
        Sam_mask  = randsample([1:numel([blastMatch_reduced.ident])],500);
        sample(i).masks = Sam_mask;
        sample(i).ident = [blastMatch_reduced(Sam_mask).ident];
        sample(i).LengthQuery = LengthQuery(Sam_mask);
        sample(i).LengthSubject = LengthSubject(Sam_mask);
        
        [testh,testp] = kstest2(sample(i).ident,[blastMatch_reduced(:).ident]);
        test(i).ks_decision = testh;
        test(i).ks_p = testp;
    end
    
    sum([test(:).ks_decision])
    
    % Plot the different sets
    
    figure(10)
    for i=1:10
        k(i)=histogram(sample(i).ident,'BinWidth',h.BinWidth);%,'FaceColor','none');
        hold on
    end
    
    figure(11)
    h2=bar(h.BinEdges(1:end-1)+h.BinWidth/2,2*h.Values/max(h.Values),'FaceColor','none','BarWidth',1); hold on
    ylim([0 2.1]);
    plot([mean([blastMatch_reduced(:).ident]) mean([blastMatch_reduced(:).ident])],[0 300],'LineWidth',3)
    
    for i=1:5
        figure(11)
        bar(k(i).BinEdges(1:end-1)+k(i).BinWidth/2,k(i).Values/max(k(i).Values),'BarWidth',1,'FaceAlpha',1); hold on
        hold on
    end
    
    

end

