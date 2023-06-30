%
%
%%%%%%%%%%%%%%%%   Cov2Ins created by Mona  %%%%%%%%%%%%%%%%%%%%%%%%
%
%
% What this script does: It finds insertions based on coverage data, where
% hybrids were mapped (hard) to the donor 
% 
%    input:
%%%%     -- List of the accGenes of the DONOR! (e.g. accBs166_2W23.txt)
%%%%            [-- No ORI crossing of accGenes implemented]
%%%%     -- Coverage of strainsOfInterest mapped hard(!) on donor

clearvars
close all

%%
%Specify your donor species
donor = "W23"; % "W23", "Bval", "Batro", "Geo"

covPath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper_BigData/Hybrids2Donors/coverage/";
replic = "Wns0" + string([3:8]);
covFiles = replic + "20";
covSuffix = "_2W23Donor_coverage.txt";

% Insertion parameters
minCovSigma = 1.5; % minCov = meanCov - minCovSigma*stdCov (default: 1.5) 
minLen = 40; % Insertions with a length of 50 bp are able to be detected! --
            % x Needs to be bigger than 20 to avoid artefacts when W23 is
            % donor (default: 40)

% Do you wanna save the insertions.mat - file?
saveInsertion = "ON"; 
saveName = ""; % Leave "" - empty if you like it just to be 'date'_insertions.mat

%Color specification
if length(covFiles)==1; cmp = [199 6 73]/255; else; cmp = jet(length(covFiles)); end

if donor == "W23"
    accPath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/acc/";
    accList = "accBs166_2W23.txt";
    donorsize = 4027680;
    artefacts = "";
    
elseif donor == "Bval"
    accPath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/acc/";
    accList = "accBs166_2Bval.txt";
    donorsize = 4286362;
    artefacts = "H:\0_PhD\Bioinformatics\Dictionaries\artefacts\B1662Bval_arteIns.txt";
    
elseif donor == "Batro"
    accPath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/acc/";
    accList = "accBs166_2Batro.txt";
    donorsize = 4168266;
    artefacts = "H:\0_PhD\Bioinformatics\Dictionaries\artefacts\B1662Batro_arteIns.txt";
    
elseif donor == "Geo"
    accPath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/acc/";
    accList = "accBs166_2Geo.txt"; 
    donorsize = 3873116;
    artefacts = "H:\0_PhD\Bioinformatics\Dictionaries\artefacts\B1662Geo_arteIns.txt";
else
    error("Please, specify your donor! -- Is it really " + donor + "?");
end

%% Load accRegions
fid = fopen(accPath + accList);
    imp = textscan(fid,'%f %f %f','delimiter','\t');
    fclose(fid);
    accStart=imp{1}; accEnde=imp{2}; accLen = imp{3};
    clear imp

accMask = false(donorsize,1);
for i = 1 : length(accStart)
    accMask(accStart(i):accEnde(i)) = 1;
end

% Load artefacts
if artefacts~=""
    fid = fopen(artefacts); imp = textscan(fid,'%f %f', 'HeaderLines', 1); fclose(fid);
    arteStart = imp{1}; arteEnd=imp{2};
    clear imp
    arteIdx = [];
    for i = 1 : length(arteStart); arteIdx = [arteIdx arteStart(i):arteEnd(i)]; end
end

%% Load coverage

for m = 1 : size(covFiles,2)
    
    fid = fopen(covPath + covFiles(m) + covSuffix);
        imp = textscan(fid,'%s %f %f','delimiter','\t');
        fclose(fid);
        rawPos = imp{2}; rawCov=imp{3};
        clear impaccLen
           
    meanCov = mean(rawCov(rawCov>0));
    stdCov = std(rawCov(rawCov>0));
    
        if donor == "Batro" % meanCov and stdCov are problematic here, because even in the core regions coverage is really low some times
            meanCov = mean(rawCov(rawCov>100));
            stdCov = std(rawCov(rawCov>100));
        elseif donor == "Geo"
            meanCov = mean(rawCov(rawCov>200));
            stdCov = 150;
        end

    minCov = round(meanCov - minCovSigma*stdCov);
    
    % Coverage in accRegions
    accCov = rawCov(accMask);
    accPos = rawPos(accMask);
    
    figure(1); set(gcf, "Renderer", "Painters");hold on; xlim([0 donorsize])
    accPlot = plot(rawPos(accMask~=0), accMask(accMask~=0).*450, ".", "Color", [0.4 0.4 0.4]);
    covinAcc = plot(accPos, accCov, "x", "Color", [0.4 0.4 0.4 0.2] );
    minCovPlot = plot([1 donorsize], [minCov minCov], '--', "Color", cmp(m,:));
    
    % Coverage greater than minCov
    insPos = accPos(accCov>minCov);
    insCov = accCov(accCov>minCov);
        
    % Artefact filtering
    if artefacts~=""
        arteFilt = ~ismember(insPos, arteIdx);
        insPos = insPos(arteFilt);
        insCov = insCov(arteFilt);
        clear arteFilt
    end
 
    insPlot = plot(insPos, insCov, "d", "Color", cmp(m,:));
    
    if isempty(insPos)
        continue
    end
    
    % list2table
    startPos = [];
    endPos = [];
    
    startPos(1) = insPos(1);
    j=1;
    
    for i = 1 : length(insPos)-1
        if insPos(i+1)~=insPos(i)+1
            endPos(j) = insPos(i);
            startPos(j+1) = insPos(i+1);
            j = j+1;
        end
    end
    endPos(j) = insPos(end);
    
    
    % Length filtering
    insLen = endPos - startPos + 1;
    lenFilt = insLen > minLen;
    
    startPos = startPos(lenFilt);
    endPos = endPos(lenFilt);
    insLen = insLen(lenFilt);
    
    % Write into struct
    if ~exist("insertions", "var")
        insertions = struct("Start", num2cell(startPos), "Ende", num2cell(endPos), "Length", num2cell(insLen), ...
            "minLen", num2cell(minLen), "minCov", num2cell(minCov), "Sample", covFiles(m));
    else
        insAdd = struct("Start", num2cell(startPos), "Ende", num2cell(endPos), "Length", num2cell(insLen), ...
            "minLen", num2cell(minLen), "minCov", num2cell(minCov), "Sample", covFiles(m));
        insertions = [insertions insAdd];
        clear insAdd
    end
    
end
figure(1);
legend([accPlot covinAcc minCovPlot insPlot], ["Accessory Genome", "Coverage in acc regions", "minCov - Line", "Potential insertions"]);

if ~exist("insertions", "var"); insertions = struct("Start", [], "Ende", [], "Length",[],"minLen",[], "minCov", [], "Sample", []); end

if saveInsertion == "ON" && exist("insertions", "var")
    if saveName == ""; saveName = datestr(now, "yyyymmdd")+ "_insertions"; end
    save(covPath + saveName , "insertions");
end
