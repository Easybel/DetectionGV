function [blastMatch,blastMatch_reduced,blastMatch_mm] = B_SortBlastHits(blastList,BedQuery,BedSubject,chrSubject,MultiMapMode,FixMode,excType_inQuery,excType_inSubject,excLfracQuery,excLfracSubject)

%% Part 1: where the blast list is sorted and all information are collected
%%%% -> no information is excluded except for excludeType_inQuery
%%%% (and of course all genes that had no hit are lost here)


fid = fopen(blastList); blastIn = textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f','delimiter','\t');
fclose(fid);
blast.locustag_in_Query=blastIn{1}; blast.ident=round(blastIn{3}/100,3); blast.startQuery=blastIn{7}; blast.endeQuery=blastIn{8}; blast.startSubject=blastIn{9};
blast.endeSubject=blastIn{10};

bedQuery = load(BedQuery);
bedQuery = bedQuery.uniqList;

if ~isempty(BedSubject)
    bedSubject = load(BedSubject);
    bedSubject = bedSubject.uniqList;
end

% here, specify the fields that you are interested in in the blastMatch
% file - here a lot of data is kept to track the quality of the blast
temp_string = strings(numel(blast.locustag_in_Query),1);
ts = num2cell(temp_string);
temp_NAN = nan(numel(blast.locustag_in_Query),1);
tn = num2cell(temp_NAN);
blastMatch = struct('locustag_in_Query',ts,'ident',tn,'startQuery',tn,'endeQuery',tn,'startSubject',tn,'endeSubject',tn,'LfracQuery',tn,...
    'otherHits',[],'otherHitsNum',[],'type_in_Query',ts,'GN_in_Query',ts,'product_in_Query',ts,...
    'LfracSubject',tn,'locustag_in_Subject',ts,'type_in_Subject',ts,'GN_in_Subject',ts,'product_in_Subject',ts);

searchFields = fieldnames(blastMatch);
basicFields = searchFields(1:6);

% blastMatch creates a unique list of all unique blast results -
% genes that didnt match do not appear and each other gene appears in one
% line
m=0;
for i=1:numel(blast.ident)
    
    % where can the gene be found in the Query bed file?
    idx_inBedQuery = find(strcmp(blast.locustag_in_Query{i},[bedQuery.locus_tag]));
    
    % unfortunetely it can happen, that the BSU name is not unique in
    % the Query bed file - take the first entry
    if numel(idx_inBedQuery)>1
        idx_inBedQuery = idx_inBedQuery(1);
        
    elseif isempty(idx_inBedQuery)
        continue
    end
    
    % skip the entry if it is of a kind that is not of interest: e.g. RNA
    if ismember(bedQuery(idx_inBedQuery).type2,excType_inQuery)
        continue
    end
    
    % the first condition always gives the first blast hit
    if i==1 || ~strcmp(blast.locustag_in_Query{i},blast.locustag_in_Query{i-1})
        m = m+1;
        
        blastMatch(m).locustag_in_Query = string(blast.locustag_in_Query(i));
        blastMatch(m).ident             = blast.ident(i);
        blastMatch(m).startQuery        = blast.startQuery(i);
        blastMatch(m).endeQuery         = blast.endeQuery(i);
        blastMatch(m).startSubject      = blast.startSubject(i);
        blastMatch(m).endeSubject       = blast.endeSubject(i);
        
        % calculate the fraction fo the Query gene that fits
        blastMatch(m).LfracQuery = round((blast.endeQuery(i)- blast.startQuery(i) +1) / (abs(bedQuery(idx_inBedQuery).ende - bedQuery(idx_inBedQuery).start)+1),3);
        
        % what is the product in Query type and Name
        extraFields  = ["type2","Name","product"];
        extraFields2Qu = ["type","GN","product"] + "_in_Query";
        
        for k=1:numel(extraFields)
            if ~isempty(bedQuery(idx_inBedQuery).(extraFields(k)))
                blastMatch(m).(extraFields2Qu(k)) = string(bedQuery(idx_inBedQuery).(extraFields(k)));
            end
        end
        
        clear c
        % this condition collects the most important information for the next best hits in an additional structure
    else
        for j=1:numel(basicFields)
            c.(basicFields{j}) = blast.(basicFields{j})(i);
        end
        % calculate the fraction fo the Query gene that fits
        c.LfracQuery = round((blast.endeQuery(i)- blast.startQuery(i) +1) / (abs(bedQuery(idx_inBedQuery).ende - bedQuery(idx_inBedQuery).start)+1),3);
        blastMatch(m).otherHits = [blastMatch(m).otherHits c];
        
    end
    blastMatch(m).otherHitsNum = size(blastMatch(m).otherHits,2);
end

% the structure was intitialized to big and is reshaped here
blastMatch = blastMatch(1:m);

%% loop over all entries in blastMatch and find the according genes in the Subject
% --> add their product and old_locus_tag
if FixMode ~= "ON"
    
    extraFields = ["locus_tag" extraFields];
    extraFields2Sub = ["locustag","type","GN","product"] + "_in_Subject";
    
    for i=1:numel(blastMatch)
        start = min(blastMatch(i).startSubject,blastMatch(i).endeSubject);
        ende = max(blastMatch(i).startSubject,blastMatch(i).endeSubject);
        
        %this covers the case where start and ende overlap a gene in S or are
        %equally sized
        idx_bigORequal = find(start <= [bedSubject(:).start] & ende >= [bedSubject(:).ende]);
        % this covers the case where start is after a gene start
        idx_small1 = find((start > [bedSubject(:).start] & start < [bedSubject(:).ende] & ende >= [bedSubject(:).ende]) | ...
            (start > [bedSubject(:).start] & ende < [bedSubject(:).ende]));
        idx_small2 = find((start >= [bedSubject(:).start] & ende < [bedSubject(:).ende]) | ...
            (start < [bedSubject(:).start] & ende < [bedSubject(:).ende] & ende > [bedSubject(:).start]));
        
        idx = unique([idx_bigORequal idx_small1 idx_small2]);
        
        for j=1:numel(idx)
            
            for k=1:numel(extraFields)
                
                if ~isempty(bedSubject(idx(j)).(extraFields(k)))
                    blastMatch(i).(extraFields2Sub(k))(j) = bedSubject(idx(j)).(extraFields(k));
                end
            end
            
            % now calculate the percentage of the length that is hit
            % --> this is done by multiplying two masks mith the genes onto
            % each other and counting the matching positions
            maskHit = zeros(chrSubject,1);
            maskHit(start:ende) = 1;
            maskSGene = zeros(chrSubject,1);
            maskSGene(bedSubject(idx(j)).start:bedSubject(idx(j)).ende) = 1;
            percLHit = round(sum(maskHit.*maskSGene) / sum(maskSGene),3);
            
            blastMatch(i).LfracSubject(j) = percLHit;
        end
        
    end
    
    
    % do the same for the otherHits
    for i=1:numel(blastMatch)
        if ~isempty(blastMatch(i).otherHits)
            for k =1:numel([blastMatch(i).otherHits(:).startSubject])
                start = min(blastMatch(i).otherHits(k).startSubject,blastMatch(i).otherHits(k).endeSubject);
                ende = max(blastMatch(i).otherHits(k).startSubject,blastMatch(i).otherHits(k).endeSubject);
                
                %this covers the case where start and ende overlap a gene in S or are
                %equally sized
                idx_bigORequal = find(start <= [bedSubject(:).start] & ende >= [bedSubject(:).ende]);
                % this covers the case where start is after a gene start
                idx_small1 = find((start > [bedSubject(:).start] & start < [bedSubject(:).ende] & ende >= [bedSubject(:).ende]) | ...
                    (start > [bedSubject(:).start] & ende < [bedSubject(:).ende]));
                idx_small2 = find((start >= [bedSubject(:).start] & ende < [bedSubject(:).ende]) | ...
                    (start < [bedSubject(:).start] & ende < [bedSubject(:).ende] & ende > [bedSubject(:).start]));
                
                idx = unique([idx_bigORequal idx_small1 idx_small2]);
                
                for j=1:numel(idx)
                    
                    for l=1:numel(extraFields)
                        
                        if ~isempty(bedSubject(idx(j)).(extraFields(l)))
                            blastMatch(i).otherHits(k).(extraFields2Sub(l))(j) = bedSubject(idx(j)).(extraFields(l));
                        end
                    end
                    
                    % now calculate the percentage of the length that is hit
                    % --> this is done by multiplying two masks mith the genes onto
                    % each other and counting the matching positions
                    maskHit = zeros(chrSubject,1);
                    maskHit(start:ende) = 1;
                    maskSGene = zeros(chrSubject,1);
                    maskSGene(bedSubject(idx(j)).start:bedSubject(idx(j)).ende) = 1;
                    percLHit = round(sum(maskHit.*maskSGene) / sum(maskSGene),3);
                    
                    blastMatch(i).otherHits(k).LfracSubject(j) = percLHit;
                end
            end
        end
    end
end
%% Make the blastMatch slimmer -- according to predefined settings


blastMatch_reduced = blastMatch;
blastMatch_mm = [];

if FixMode ~= "ON"
    % 4. reduce the entries of higher dimension to dim=1 -- some blasts hit more than one gene
    % --> look at LfracSubject and take the gene where this is maximal
    
    for i=1:numel(blastMatch)
        %         a = []; b = []; c = []; d = [];
        if numel(blastMatch_reduced(i).LfracSubject)>1
            BestHit = find(blastMatch_reduced(i).LfracSubject == max(blastMatch_reduced(i).LfracSubject));
            
            if numel(BestHit)==1
                bestH = BestHit;
            else
                bestH = BestHit(1);
            end
            
            for k=1:numel(extraFields2Sub)
                if numel(blastMatch_reduced(i).(extraFields2Sub(k))) == numel(blastMatch_reduced(i).LfracSubject)
                    blastMatch_reduced(i).(extraFields2Sub(k)) = blastMatch_reduced(i).(extraFields2Sub(k))(bestH);
                end
            end
            
            blastMatch_reduced(i).LfracSubject = blastMatch_reduced(i).LfracSubject(BestHit(1));
            
            %%%%%%% THIS IS %%%%%%% OPTIONAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if two genes are both the best match then write them both
            %for j=1:numel(BestHit)
            %    b = [b  ' | ' blastMatch_reduced(i).type_in_Subject{BestHit(j)}];
            %    c = [c  ' | ' blastMatch_reduced(i).GN_in_Subject{BestHit(j)}];
            %    d = [d  ' | ' blastMatch_reduced(i).locustag_in_Subject{BestHit(j)}];
            %end
            
            %blastMatch_reduced(i).type_in_Subject = b;
            %blastMatch_reduced(i).GN_in_Subject = c;
            %blastMatch_reduced(i).locustag_in_Subject = d;
            %%%%%%% THIS IS %%%%%%% OPTIONAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear BestHit
            
        end
    end
    
    % 2. kick out those where the length does not match enough
    keep_Length = ([blastMatch_reduced(:).LfracQuery] >= excLfracQuery & [blastMatch_reduced(:).LfracSubject] >= excLfracSubject  );
    blastMatch_reduced = blastMatch_reduced(keep_Length);
    
    if MultiMapMode == "ON"
        
        % 1b. throw out everything where there is only 1 hit and it is the same as
        % the query gene
        query_is_subject = [blastMatch_reduced(:).locustag_in_Query]==[blastMatch_reduced(:).locustag_in_Subject];
        keep_mm          = [blastMatch_reduced(:).otherHitsNum]>1;
        
        blastMatch_mm = blastMatch_reduced(keep_mm);
        blastMatch_mm = rmfield(blastMatch_mm,{'ident','startSubject','endeSubject','LfracSubject','locustag_in_Subject','type_in_Subject','product_in_Subject','GN_in_Subject'});
        % kick out those, where no gene is found in the otherHits
        
        
        for i=1:numel([blastMatch_mm(:).startQuery])
            
            if numel(fieldnames(blastMatch_mm(i).otherHits))< 8
                blastMatch_mm(i).otherHits          = [];
                blastMatch_mm(i).otherHitsNum       = 0;
                blastMatch_mm(i).maxHitLfracQuery  = 0;
            else
                hasAName   = cellfun(@(x) ~isempty(x),{blastMatch_mm(i).otherHits(:).locustag_in_Subject});
                LfracCheck = [blastMatch_mm(i).otherHits(:).LfracQuery] >= excLfracQuery;
                
                blastMatch_mm(i).otherHits          = blastMatch_mm(i).otherHits(hasAName & LfracCheck);
                blastMatch_mm(i).otherHitsNum       = sum(hasAName & LfracCheck);
                if ~isempty(max([blastMatch_mm(i).otherHits(:).LfracQuery]))
                    blastMatch_mm(i).maxHitLfracQuery   = max([blastMatch_mm(i).otherHits(:).LfracQuery]);
                else
                    blastMatch_mm(i).maxHitLfracQuery   = 0;
                end
            end
        end
        
        blastMatch_mm = blastMatch_mm([blastMatch_mm.maxHitLfracQuery] > 0);
        % 2. kick out those where the length does not match enough
        keep_Length = ([blastMatch_reduced(:).LfracQuery] >= excLfracQuery & [blastMatch_reduced(:).LfracSubject] >= excLfracSubject);
        blastMatch_reduced = blastMatch_reduced(keep_Length);
        
        
    end
    
    
    % 1. Throw out everything of type (excType_inSubject)
    
    keep_TypeSubject = cellfun(@(x) ~strcmp(excType_inSubject,x),{blastMatch_reduced(:).type_in_Subject});
    blastMatch_reduced = blastMatch_reduced(keep_TypeSubject);
    
end
% 3. remove other Hits column, as it isn't interesting here, where the
% second. third ... best hit of a gene would have been
blastMatch_reduced = rmfield(blastMatch_reduced,'otherHits');

end
