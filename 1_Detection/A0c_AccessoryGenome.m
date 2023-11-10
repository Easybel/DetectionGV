% % ------ % % ------ New, accgenome.m -@author:MF ------ % % ------ % %
% % - This script needs a long time, since it's working - % % ------ % %
% % - with the *_bcf.vcf - file (~4M rows * 10 columns) - % % ------ % %
clearvars

% I only need the *bcf.vcf - file, and some information about the parameter
% than we can start:

    path = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper_BigData/Donor2Bs166/"; 
    fileName   = "Bmoj_2NCe";
    fileSuffix = "_bcf.vcf";

    % genome size of the .dict that .vcf was mapped to = number of letters in .fasta file
    dictsize = 4215607; % B.sub Bs166
    %dictsize = 4027680; % B.spiz W23 WT
    
    maxDP = 50;     % Everything covered >= maxDP, is not included in the 
                    % accessory genome
    minLength = 0;    % Only regions longer minLength bps are considered
                    
    saveName = "acc" + fileName;
    savePath = path;
    
%
% % That's all, thank you.
%

% Read in bcf.vcf - file
fid = fopen(path + fileName + fileSuffix);
imp = textscan(fid,'%*s %d . %*s %*s %*f . %s %*s %*s', 'CommentStyle','#'); % Columns with %* will not be read in
fclose(fid);
pos = imp{1};
info = imp{2};
clear imp

% Filter INDELS
filterIndels = ~contains(info, "INDEL");
posAll = pos(filterIndels);
infoAll = info(filterIndels);

% Extract DP info
infosplit = cellfun(@(s) split(s, ["DP=" ";I16="]), infoAll, 'UniformOutput', false);
DPcell = cellfun(@(s) str2num(s{2,1}), infosplit, 'UniformOutput', false);
DP = cell2mat(DPcell);


% 1. Find positions with DP < maxDP

DPmask = DP < maxDP;
posDP = posAll(DPmask);

% 2. Find positions with no entry in the bcf - no read mapped, position
% with no cover --> posMissingMask == 1

posRecip = [1 : dictsize];
posMissingMask = ~ismember(posRecip, posAll);

% include 1.: set all pos where DP < 50: posMissingMask == 1
posMissingMask(posDP) = 1;

% Translate into positions
posMissing = posRecip(posMissingMask);

% List 2 Table
startIdx = find(posMissing - [1 posMissing(1:end - 1)] ~= 1);
endIdx = [startIdx(2:end) - 1, numel(posMissing)];

start = posMissing(startIdx);
fin = posMissing(endIdx);
totalLength = (fin - start + 1); 

% Length Filtering
lengthFilt = totalLength > minLength;

startFilt = start(lengthFilt);
endFilt = fin(lengthFilt);
totalLengthFilt = totalLength(lengthFilt);

% Saving the filtered list
accList = [startFilt; endFilt; totalLengthFilt];
accOutput = fopen(savePath + saveName + "_allBP.txt", "w");
fprintf(accOutput, '%d %d %d\n', accList);
fclose(accOutput);

save(savePath + saveName + "_allBP.mat","accList")

% 
accessoryFraction = sum(totalLengthFilt)/dictsize;
disp("For " + fileName + ", "+ num2str(accessoryFraction*100, "%0.1f")+ " % of the genome is accessory.");
