%
%
% I want all your INDELS !! and check if they are on the masterlist
% 
clearvars; close all

path = "H:\fitnessDistribution\2022_withSelection\DNASeq\mappedData\";
path = "H:\phenotypeProject\DNAseq\";
buzzword = "_";
qualFilt = 0; % 0 if you do not want to filter for phredQual
donor = "Bval"; 

    % Import artefact-indels
    arteIndels = "H:\0_PhD\Bioinformatics\Dictionaries\artefacts\Artefacts_Bs1662BsNCe.vcf";
    fid = fopen(arteIndels);
    imp = textscan(fid, '%s %f %s %s %s %f %s %s %s %s', 'CommentStyle', '#');
    impPos = imp{2}; impInfo = imp{8}; 
    arteIndelPos = impPos(contains(impInfo, "INDEL"));
    fclose(fid);
    clear imp*;

    if donor == "Bval"
        mlBcfc = "H:\masterlistMapping\mlBval2Bs166NCe\Bval2Bs166NCe_bcfcall.vcf";
    end
    
    fid = fopen(mlBcfc);
    imp = textscan(fid,'%*s %d . %s %s %f . %s %*s %s', 'CommentStyle','#'); % Columns with %* will not be read in
    fclose(fid);
    impPos = imp{1}; impInfo = imp{5}; 
    mlIndelPos = impPos(contains(impInfo, "INDEL"));
    
allFiles = dir(path);
fileNames = {allFiles([allFiles(:).isdir]==0 & contains({allFiles(:).name}, "bcfcall")& contains({allFiles(:).name}, buzzword)).name};

indelSummary = [];
for m = 1 : size(fileNames,2)
    
        fid = fopen(path+fileNames{m});
        imp = textscan(fid,'%*s %d . %s %s %f . %s %*s %s', 'CommentStyle','#'); % Columns with %* will not be read in
        fclose(fid);

        % In the struct snp, you will now write Position, Reference, Alternate, 
        % Info and Counts :
        snp = cell2struct(horzcat(num2cell(imp{1}),imp{2:3}, num2cell(imp{4}), imp{5:end}), ...
            ["pos","ref","alt","qual","info","counts"],2);
        clear imp fid;

        % Find indels
        indelFind = contains({snp.info}, "INDEL");
        indelsAll = snp(indelFind);
        
        % Filter qual + exclude artefacts
        indelsFilt{m} = indelsAll([indelsAll.qual]>qualFilt & ~ismember([indelsAll.pos], arteIndelPos));
        
        % Check if they are on the masterlist (only position)
        potentialMLindel = num2cell(ismember([indelsFilt{m}.pos], mlIndelPos)');
        sampleName = repmat(cellstr(fileNames{m}),size(indelsFilt{m},1),1);
        
        [indelsFilt{m}.onML]=potentialMLindel{:};
        [indelsFilt{m}.sample]=sampleName{:};
        
        indelSummary = [indelSummary; indelsFilt{m}];
end