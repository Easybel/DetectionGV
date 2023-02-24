 % %   %   %   %  % %                            % %  %   %   %   % %
% %   %   %   %  % %   Master List Filtering . m  % %  %   %   %   % %
 % %   %   %   %  % %                            % %  %   %   %   % %
clearvars

     %                                              %
    % %   Help me doing the right thing for you:   % %
     %                                              %
     
%%% Where are your master list mapping results?
inppath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper_BigData/Donor2Bs166/"; 
inpname = "Bmoj_2NCe"; 

inpsuffix = "_bcfcall.vcf";

%%% Where do you want to save the master list?
savepath = inppath;
savename = "mlBmoj2Bs166NCe_v1";

% What filtering parameters would you like to set ?
minTotalCount = 50;     % What is your total count cut off ? (default: 50)
totalCountsMode = "DP"; % Do you like the raw read depth DP or the sum over 
                        % all allelic depths AD be taken as your total 
                        % counts ? ("DP" / "AD")
minQual = 50;           % What is the minimum phred quality ? (default: 50)
unambiFreq = 90;        % [in %] When is a variant unambigious ? (default : 90)
excludeIndels = "ON";   % Do you want to keep indels?

                 %   Thanks for your help!
                % %  
                 %   From here, I can be on my own. 

inpu = inppath + inpname + inpsuffix;

mlSummary = struct("masterlist", [], ...
        "inputfile",cellstr(inpname + inpsuffix), "outputfile", cellstr(savename+".txt"),...
        "filterDate", datestr(now, 'dd/mm/yyyy'), "totalCountMode", totalCountsMode, ...
        "minTotalCounts", minTotalCount, "minQual", minQual,...
        "unambigiousFreq", unambiFreq, "LAfiltered", false, ...
        "unfiltered", [] ... 
        );
%% Read data from bcf-call
    fprintf("Read bcf-call '%s'\n", inpu);

    fid = fopen(inpu);
    imp = textscan(fid,'%*s %d . %s %s %f . %s %*s %s', 'CommentStyle','#'); % Columns with %* will not be read in
    fclose(fid);

    % In the struct snp, you will now write Position, Reference, Alternate, 
    % Info and Counts :
    mlsnp = cell2struct(horzcat(num2cell(imp{1}),imp{2:3}, num2cell(imp{4}), imp{5:end}), ...
        ["pos","ref","alt","qual","info","counts"],2);
    clear imp fid;

%% Convert bcf-call data to a useable data structure
    fprintf("\tConvert data...\n");

    splittedCounts = split({mlsnp.counts}',":"); %just splitting last bcfcall-column at :
    altsArr = {mlsnp.alt};    % all alternates
    altCounts = count({mlsnp.alt}',",") + 1;  % count number of alternates

    % Preallocate memory for countArrs and altArrs
    countArrs = cell(size(mlsnp));
    altArrs = cell(size(mlsnp));
    % Splitting splittedCounts and altsArr:
    for alts = 1:3
        maskPos = altCounts == alts;
        countArrs(maskPos) = num2cell(str2double(split(splittedCounts(maskPos,end),",",2)),2);
        altArrs(maskPos) = mat2cell(split(altsArr(maskPos)',",",2),ones(sum(maskPos),1));
    end
    clear splittedCounts altsArr altCounts alts mask;

    totalCounts = cellfun(@(x) sum(x), countArrs);

    % Now, decide for the alt with the higher frequency (more counts)
    altAllele = cellfun(@(x) outputs2array_LH(@()max(x(2:end)),[1 2]), countArrs, 'UniformOutput',false);
    altAllele = vertcat(altAllele{:});
    altAlleleChar = cellfun(@(alt,idx) alt{idx}, altArrs, num2cell(altAllele(:,2)), 'UniformOutput', false);
    altFrequency = altAllele(:,1) ./ totalCounts;

    clear altArrs altAllele;

    if strcmp(totalCountsMode, 'DP')
        % Extract DP info
        cellY = cellfun(@(s) strsplit(s, 'DP='), {mlsnp.info}', 'UniformOutput', false);
        cellZ = cellfun(@(s) strsplit(s{2}, ';'), cellY, 'UniformOutput', false);
        cellDP = cellfun(@(s) str2double(s{1}), cellZ, 'UniformOutput', false);
        DP = cell2mat(cellDP);
        clear cell*
    end


%% Creating filters ...
    fprintf("\tCreate filter...\n");
    
    filterIndels = ~contains({mlsnp.info},"INDEL");
    filterQual = [mlsnp.qual]' >= minQual; 
    filterUnambi = altFrequency >= unambiFreq / 100; 

    if strcmp(totalCountsMode, 'AD')
        filterTotalCounts = totalCounts >= minTotalCount;
    elseif strcmp(totalCountsMode, 'DP')
        filterTotalCounts = DP >= minTotalCount;
    else
        error("Your total count filtering fails! Please set your totalCountsMode to AD or DP.");
    end

    % Merging all filters ...
    if excludeIndels == "ON"
    filterAll = filterIndels' & filterTotalCounts & filterQual & filterUnambi;
    else
    filterAll = filterTotalCounts & filterQual & filterUnambi;
    end


%% Write in mlSummary 
    fprintf("\tWrite to summary...\n");
    
    % Write unfiltered data in Summary
    mlSummary.unfiltered = ...
        struct("pos",{mlsnp.pos}', "DP", num2cell(DP), "ref",{mlsnp.ref}', "alt",{mlsnp.alt}', ...
        "varcounts",countArrs,"varfreq", num2cell(altFrequency), "QUAL", {mlsnp.qual}', "filterIndels",num2cell(filterIndels)',"filterQual", num2cell(filterQual), ...
        "filterTotalCounts",num2cell(filterTotalCounts), "filterUnambi",num2cell(filterUnambi), ...
        "filterAll",num2cell(filterAll)); 

    % Filter data
    mlsnp_filt = mlsnp(filterAll);

    % Create masterList in mlSummary
    mlSummary.masterlist = vertcat([mlsnp_filt.pos], string({mlsnp_filt.ref}), string(altAlleleChar(filterAll))')';

    clear countArrs altAlleleChar filterIndels filterTotalCounts ...
        filterCoverage filterAll snp;
    
        
%% Save data 

    disp("Save master list (*.txt) ...");

    fileID = fopen(savepath + savename + ".txt",'w');
    fprintf(fileID,"%s %s %s\n", [mlSummary.masterlist]' );
    fclose(fileID);

    fprintf("Save FilterSummary and MasterList to '%s.mat'...\n", savepath + savename);
    save(savepath + savename, 'mlSummary');

    disp("Done.");