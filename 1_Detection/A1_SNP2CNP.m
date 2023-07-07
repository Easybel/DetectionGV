          % %                                       % %
         % % %                                     % % % 
        % % % %    SNP2CNP created by MonaIsa     % % % %
         % % %                                     % % %
          % %                                       % % 

% This script takes the filtered variants and finds the replacements by 
% detecting clusters of donor alleles. 

% input: 
% -- SNPSummary/ MutSummary (which is output from A0_VariantFiltering.m)
% -- masterlist containing SNPwise differences between recipient and donor
%    (obtained with A0b_MasterListFiltering.m)
% -- accessory genome of recipient with respect to donor 
%    (obtained with A0c_AccessoryGenome.m)
% -- list with multimapper regions (obtained with A0d_Multimapper.m)
% -- list of SNP artefacts that are found with the A0_VariantFiltering.m
%    script together with the recipient lab strains aligned to its
%    own reference

% output: CNPSummarya .mat
% -->>> structure containing all detected replacements (Adist, Mdist, Cdist); 
%       de novo mutations, de novo indels; etc.
% C ->>> this field contains cluster start and end positions in C(1,:), C(2,:)
% --->>> number of total masterlist SNPs in cluster region C(3,:)
% --->>> number of undetected mlSNPs in cluster C(4,:)
%
% Calculating cluster lengths (distances between cluster start and end)
% Cdist: minimal length, between first and last detected mlSNP
% Mdist: maximal length, including the identical regions befor and after  
% ------- from next mlSNP that is undetected before and after cluster
% Adist: average(a))
% ------- average between last undetected and first detected SNP in
%         cluster, at start and end
%
% identity, mlSNPs in cluster with Adist
% ident{m} = 1 - (# mlSNPs) / Adist
%
% transferred core genome per sample: evaluated with Adist

clearvars
%%
                    %                      %
                   % %    Define input    % %
                    %                      %

% select the used donor species for masterlist and accessory genome
donor = "Bspiz"; % "Bval", "Batro", "Geo", "Bmoj"

pathLists = "/DetectionGV/dictionaries_Bacillus/";
if strcmp(donor, "Bval")
    % Bval donor
    masterlist = pathLists + "ml/" + "mlBval2Bs166NCe_v1.txt";
    accgenome  = pathLists + "acc/" + "accBval_2NCe.txt";
elseif strcmp(donor, "Bspiz")
    % Bspiz donor
    masterlist = pathLists + "ml/" + "mlW232Bs166NCe_v1.txt";
    accgenome  = pathLists + "acc/" +"accW23_2NCe.txt";
elseif strcmp(donor, "Batro")
    % BA donor
    masterlist = pathLists + "ml/" + "mlBatro2Bs166NCe_v1.txt";
    accgenome  = pathLists + "acc/" +"accBatro_2NCe.txt";
elseif strcmp(donor, "Geo")
    % Geo donor
    masterlist = pathLists + "ml/" + "mlGeo2Bs166NCe_v1.txt";
    accgenome  = pathLists + "acc/" +"accGeo_2NCe.txt";
elseif strcmp(donor, "Bmoj")
    % Bmoj donor
    masterlist = pathLists + "ml/" + "mlBmoj2Bs166NCe_v1.txt";
    accgenome  = pathLists + "acc/" +"accBmoj_2NCe.txt";
else
    error("Something went wrong with your donor declaration.. Is your donor really %s?", donor);
end

%define the multimapper list and the list of SNP artefacts
mmlist    = pathLists + "Bs166NCe_mm.txt";
artefacts = pathLists + "Bs166NCe_ArteSNPs.vcf";

recipsize = 4215607;

%
% you have to pick one type and cant mix them!
SNPSource = "MutSummary";     % Where do your variant lists come from? 
                              % "SNPSummary" (default)
                              % "MutSummary" (includes Indels) 
                              % "IndvMutLists" (old: is a .txt file)
                             
% You just need to specify samplenames if you use "IndvMutLists"
% otherwise leave empty to choose ALL:
samplenames = [];

% Path of the SNPSummary/ies or IndvMutLists of the evolved strains:
SNPPath = "../IN";

% Turn "ON" to save the CNPSummary !
saveCNPSummary = "ON"; 

if any(strcmp(SNPSource, ["SNPSummary", "MutSummary"]))
    % You can give more than {1} SNPSummary as input
    SNPName{1} = "XXX.mat";


elseif strcmp(SNPSource, "IndvMutLists")
    % Default setting: The IndvMutLists are saved as the samplenames 
    % plus "_IndvMutList.txt"
    SNPName = samplenames + "_IndvMutList.txt";
    
else
    error("Please check your SNPSource! It needs to be 'SNPSummary', 'MutSummary' or 'IndvMutLists'.");
end

% Where do you want to save the output -- CNPSummary ?
savepath = SNPPath;

% Define cluster finding parameters
cms = 5;    % cms - 1: number of allowed missing snps!!
            % for cms and more, the cluster ends

%%
  %                                                                  %
 % %                                                                % %
% % %   Now, just lean back and let me do the rest for you ...     % % %
 % %                                                                % %
  %                                                                  %
  
% Load variant lists into allSNPs
allSNPs = [];
if any(strcmp(SNPSource, ["SNPSummary", "MutSummary"]))
    
    for i = 1 : length(SNPName)
        Data = load(SNPPath + SNPName{i});
        allSNPs = [allSNPs Data.SNPSummary];
    end
    if isempty(samplenames)
        samplenames = [allSNPs.Sample];
    end
    
elseif strcmp(SNPSource, "IndvMutLists")
    
    for i = 1 : length(SNPName)
        fid = fopen(SNPPath + SNPName(i));
            imp = textscan(fid, '%f %s %s');
            allSNPs(i).Sample = samplenames(i);
            allSNPs(i).IndvMutList = [imp{1}, string(imp{2}), string(imp{3})];
            clear imp;
        fclose(fid);
    end
    disp("You are using IndvMutList.txt - files to analyse your data.");
    disp("These are obsolete (old filtering)!! Check out the new SNPSummaries!:)");
    
end
  
% Preallocating your variables
indel  = cell(numel(samplenames),1);     % this is only needed for MutSummary that contains indels
denovo = cell(numel(samplenames),1);     % SNPs that do not fit the masterlist and are no artefacts
C      = cell(numel(samplenames),1);     % detected clusters
SPI    = cell(numel(samplenames),1);     % SNPs that fit the masterlist but are not in cluster (You need at least two mlSNPs to have a cluster)
Cdist  = cell(length(C),1);              % minimum length of the detected cluster (from first donor allele SNP to the last)
Mdist  = cell(length(C),1);              % maximum length (including the distance to the next mlsnp or to the start of the next accessory part)
Adist  = cell(length(C),1);              % average length ((maxlen-minlen)/2+minlen)
ident  = cell(length(C),1);
transfer = zeros(1,length(C));
ORI_CROSSING = zeros(length(C),1);

% Load recipient/donor specific master list and write position in ml.pos and
% the alternate base in ml.atcg (ml is a struct)
fid = fopen(masterlist);
imp = textscan(fid,'%f %s %s');
fclose(fid);
ml.pos = imp{1};
ml.atcg = imp{3};
clear imp

% Load recipient/donor specific accessory genome in acc (as a structure)
fid = fopen(accgenome);
imp = textscan(fid, '%f %f %f', 'delimiter', '\t');
acc.start = imp{1};
acc.edge = imp{2};
fclose(fid);
clear imp

% Load the multimapper list 
fid = fopen(mmlist);
imp = textscan(fid, '%f %f %f');
mm.start = imp{1};
mm.edge = imp{2};
clear imp;

% Merge with acc.list for later use:
exclRegion.start = vertcat(acc.start, mm.start);
exclRegion.edge = vertcat(acc.edge, mm.edge);

% Create excl_mask that will exclude all acc. and mm regions
excl_mask = [];
for i = 1 : length(exclRegion.start)
    excl_mask = [excl_mask exclRegion.start(i):exclRegion.edge(i)];
end 
excl_mask = unique(sort(excl_mask));

% Read in artefacts.vcf (recipient reads aligned to its own dictionary)
fid = fopen(artefacts);
imp = textscan(fid, '%s %f %s %s %s %f %s %s %s %s', 'HeaderLines', 33);
artefact.pos = imp{2};
artefact.alt = imp{5};
fclose(fid);
clear imp; 

%% This loop goes through all samples in samplenames

for m = 1:numel(samplenames)
    
% import SNP list of the current sample
clear snp
findSample = strcmp(samplenames(m), [allSNPs.Sample]);

if length(nonzeros(findSample)) == 0
    error("Your input data does not contain sample '%s'!", samplenames(m));
elseif length(nonzeros(findSample)) > 1
    fprintf("Be aware that your sample '%s' can be found more than once in the loaded data!\n",  samplenames(m));
    fprintf("Only the first entry will be considered in the following analysis!\n");
    idx = find(findSample == 1);
    findSample = min(idx);
end

snp.pos = str2double(allSNPs(findSample).IndvMutList(:,1));
snp.ref = allSNPs(findSample).IndvMutList(:,2);
snp.atcg = allSNPs(findSample).IndvMutList(:,3);

% remove the artefacts 
for i = 1 : length(artefact.pos)
    idx = find(snp.pos == artefact.pos(i));
    if ~isempty(idx) && strcmp(artefact.alt(i), snp.atcg(idx))
        snp.pos(idx) = []; snp.atcg(idx) = []; snp.ref(idx) = [];
    end   
end 
clear idx 

if isempty(snp.pos)
    continue
end

% 
% find matching snps and masterlist snps!!
%
% first exclude all the indels and write them into indel list for later use
% (we don't want to search clusters based on indels)

indel_mask        = cellfun(@(x) numel(x),snp.ref)>1 | cellfun(@(x) numel(x),snp.atcg)>1;
if any(indel_mask)
    indel{m}.pos(1,:) = snp.pos(indel_mask)';
    indel{m}.pos(2,:) = 0;
    indel{m}.var      = snp.atcg(indel_mask)';
    indel{m}.ref      = snp.ref(indel_mask)';

    % remove indels from snp
    snp.atcg = snp.atcg(~indel_mask);
    snp.pos = snp.pos(~indel_mask);
    snp.ref = snp.ref(~indel_mask);
end

if isempty(snp.pos)
    continue
end

% % % % % % % % % % 
% % % % % % % % % % 

% Match function (MATCH_MF): which positions are on both lists? 
% -- match_pos is a mask (logicals) for the masterlist
% -- match_posidx gives the index on the snplist of the replicate

[match_pos, match_posidx] = ismember(ml.pos, snp.pos);    

% Are the detected alternate alleles the same as the expected?
match_nt = strcmp(ml.atcg(match_pos), snp.atcg(nonzeros(match_posidx)));    

paren = @(x, varargin) x(varargin{:});      % definition of an anonymous function
match = paren(ml.pos(match_pos), match_nt); % Match contains the positions on the reference genome where mlSNP was detected

match_mask = ismember(ml.pos, match);       % creates a masterlist mask with detected mlSNPs

% % % % % % % % % % 
% % % % % % % % % % 

%
% Create list of denovo mutations 
%
denovo{m}.pos(1,:) = snp.pos(~ismember(snp.pos, match))';
denovo{m}.pos(2,:) = 0;
denovo{m}.ref = snp.ref(~ismember(snp.pos, match))';
denovo{m}.var = snp.atcg(~ismember(snp.pos, match))';
%
    
if isempty(match)
    continue
end
   
%
% Finding Cluster
%  
% Finding the start of the clusters:
% For this, the masterlist mask ($match_mask) is shifted $cms-times to the right and
% added to the non shifted mask, $sta, every time. For every start of a cluster (a detected
% masterlist SNP after $cms of undetected SNPs), $stat contains a 0 followed
% by a 1. Then, $sta is divided by $sta shifted by 1 to the right ($sta2)
% sucht that every start position with respect to the masterlist is INF. 
% All other positions are an INT or NaN

clear sta sta2 start
sta = match_mask(1:end)'; 
for i = 1:cms-1
    sta = [match_mask(end-i+1:end)' match_mask(1:end-i)'] + sta; 
end
sta2 = sta(1:end)./[sta(end:end) sta(1:end-1)];
start = ml.pos(isinf(sta2));

% Finding the ends of the clusters(edge):
% The end of the cluster is found the same way as the start, just that the 
% mask is shifted to the left and summed up ($ed).Then, $ed is divided by $
% ed2 and the end positions of clusters are saved in $edge

clear ed ed2 edge
ed = match_mask(1:end)';
for i = 1:cms-1
    ed = [match_mask(1+i:end)' match_mask(1:i)'] + ed; 
end 
ed2 = ed(1:end) ./ [ed(2:end) ed(1:1)];
edge = ml.pos(isinf(ed2)); %gives only the ml positions where ed2 equal inf

% Tidy up
clear sta sta2 ed ed2

%
% Test, if a cluster is going over the ORI , if yes , split first
%

if edge(1) < start(1)
    start = vertcat(0, start);
    edge = vertcat(edge, recipsize);
    ORI_CROSSING(m) = 1;
    disp(['Your replicate ', samplenames{m}, ' has an ORI crossing cluster, dude!']);
end

% All clusters with length 1 are counted as SPIs
% -- longer clusters are kept 
SPI_mask = ismember(start, edge);
SPI{m,1} = start(SPI_mask);
start = start(~SPI_mask);
edge = edge(~SPI_mask);
clear SPI_mask

if isempty(start) && isempty(edge)
    continue
end

% Create a clustermask with all positions within a cluster now, you will
% need it for the exclusion of the accessory parts and to find if denovos
% are within clusters or outside
clusterPos = [];
for i = 1 : length(start)
    clusterPos = [clusterPos start(i):edge(i)];
end 

% From the clusters, multimapping regions and accessory parts are excluded
% and cut out here:
excl_filter = ~ismember(clusterPos, excl_mask);
clusterPos_filt = clusterPos(excl_filter);
clear clusterPos acc_filt

% Are denovos in a cluster or outside ?
if ~isempty(denovo{m})
denovo{m}.pos(2,:) = ismember(denovo{m}.pos(1,:), clusterPos_filt);
end
% Are indels in a cluster or outside ?
if ~isempty(indel{m})
indel{m}.pos(2,:) = ismember(indel{m}.pos(1,:), clusterPos_filt);
end

% Now go back from the list of all positions in a cluster to start/end
% position of the detected (acc.gen. and mm-regions - cleaned) cluster
k = 1;
start_filt(k) = clusterPos_filt(1);

for i = 1 : length(clusterPos_filt) - 1
   if clusterPos_filt(i+1) ~= clusterPos_filt(i)+1
       edge_filt(k) = clusterPos_filt(i);
       start_filt(k+1) = clusterPos_filt(i+1);
       k = k + 1;
   end 
end
edge_filt(k) = clusterPos_filt(end);

start = start_filt;
edge = edge_filt;
clear start_filt edge_filt

% Starting to fill the cell C with information about the detected clusters 
C{m}(1,:) = start';
C{m}(2,:) = edge';

% Count masterlist SNPs in clusters: 
% -- all SNPs from cluster start - SNPs after cluster end
for i = 1 : length(start)

    % Number of mlSNPs in cluster
    snpno(i) = length(ml.pos(ml.pos >= start(i))) - length(ml.pos(ml.pos > edge(i)));
    C{m}(3,i) = snpno(i);

    % Number of undetected mlSNPs in cluster:
    % total mlSNPs - detected mlSNPs in region
    snpno_det(i) = length(match(match >= start(i))) - length(match(match > edge(i)));
    C{m}(4,i) = C{m}(3,i) - snpno_det(i); % number of missing snps
    clear snpno snpno_det
    
    % Calculating cluster lengths (distances between cluster start and end)
    % Cdist: minimal length, between first and last detected mlSNP
    % Mdist: maximal length, including the identical regions befor and after  
    % ------- from next mlSNP that is undetected before and after cluster
    % Adist: average(a))
    % ------- average between last undetected and first detected SNP in
    %         cluster, at start and end

    cdist = edge(i) - start(i) + 1;
    Cdist{m}(i,1) = cdist;

    % Measuring cluster lengths: To calculate the maximum distance, the next mlSNP that is not in the
    % cluster needs to be found or the next accessory part (which one of both is closer)
    mdist_back = min(ml.pos(ml.pos > edge(i)));         % find the next mlsnp behind the cluster
    if isempty(mdist_back); mdist_back = recipsize; end % if there is no, the maximum is the end of the genome
    
    mdist_front = max(ml.pos(ml.pos < start(i)));   % find the next mlsnp in front of the cluster
    if isempty(mdist_front); mdist_front = 0; end   % if there is no, the limit is the start of the genome
    
    idx_back = find(edge(i) < exclRegion.start(:,1) & mdist_back > exclRegion.start(:,1)); % see if there is an acc part between the two snps in the back
        if ~isempty(idx_back) 
            mdist_back = exclRegion.start(min(idx_back));  % if there is an acc part, set the range of mdist to the start of the accessory part
        end    
    idx_front = find(start(i) > exclRegion.edge(:,1) & mdist_front < exclRegion.edge(:,1)); % see if there is an acc part between the two snps in the front
        if ~isempty(idx_front) 
            mdist_front = exclRegion.edge(max(idx_front)); % if there is an acc part, set the range of mdist to the end of the accessory part
        end  
    mdist = (mdist_back - 1) - (mdist_front + 1) + 1;
    Mdist{m}(i,1) = mdist;
    Mdist{m}(i,2) = mdist_front + 1;
    Mdist{m}(i,3) = mdist_back - 1;
        
    % Average distance: adist = cdist + (mdist - cdist)/2;
    Adist{m}(i,2) = ceil(mdist_front + (start(i) - mdist_front)/2);
    Adist{m}(i,3) = ceil(mdist_back - (mdist_back - edge(i))/2);
    Adist{m}(i,1) = Adist{m}(i,3) - Adist{m}(i,2) + 1;

    clear mdist cdist mdist_back mdist_front
end 
   
% Now merge, the ORI crossing cluster 
if ORI_CROSSING(m)~=0
    C{m}(1,1) = C{m}(1,end); % 
    C{m}(3,1) = C{m}(3,1) + C{m}(3,end);
    C{m}(4,1) = C{m}(4,1) + C{m}(4,end);
    C{m}(:,end) = [];
    
    Adist{m}(1,1) = Adist{m}(1,1) + Adist{m}(end,1);
    Adist{m}(1,2) = Adist{m}(end,2);
    Adist{m}(end, :) = [];
    
    Mdist{m}(1,1) = Mdist{m}(1,1) + Mdist{m}(end,1);
    Mdist{m}(1,2) = Mdist{m}(end,2);
    Mdist{m}(end, :) = [];
    
    Cdist{m}(1,1) = Cdist{m}(1,1) + Cdist{m}(end,1);
    Cdist{m}(end, :) = [];
end 
    
% Calculate the identity:
% Identity = 1 - (# mlSNPs) / Adist
ident{m} = 1 - C{m}(3,:)' ./ Adist{m}(:,1);
 
% Calculate the percentage of genome replaced for Adist 
transfer(m) = sum(Adist{m}(:,1)) / recipsize * 100 ;

% Tidy up
clear i idx_back idx_front fid match_posidx match_pos match_nt
end

if ~isempty(nonzeros(transfer))
    disp(['On average, (', num2str(mean(transfer(transfer>0)), '%.2f'),' ',char(177),' ',num2str(std(transfer(transfer>0)), '%.2f'),') % of the genome is replaced (replicates with no clusters excluded).']);
end

%save 
samplenamesC = cellfun(@(x) {string(x)},cellstr(samplenames));
CNPSummary = struct('C',C,'Cdist',Cdist,'Adist',Adist,'Mdist',Mdist,'denovo',denovo,'SPI',SPI, 'Ident', ident,'ORI_Crossing',num2cell(ORI_CROSSING),'Transfer', num2cell(transfer'),'Sample', samplenamesC');
if SNPSource == "MutSummary"
    CNPSummary = struct('C',C,'Cdist',Cdist,'Adist',Adist,'Mdist',Mdist,'denovo',denovo,'indel',indel,'SPI',SPI, 'Ident', ident,'ORI_Crossing',num2cell(ORI_CROSSING),'Transfer', num2cell(transfer'),'Sample', samplenamesC');
end

if saveCNPSummary == "ON" 
    save(savepath + datestr(now, 'yyyymmdd') + "_CNPSummary.mat", 'CNPSummary')
end

% Tidy up; comment out if you like to have a closer look at variables
clear excl_mask clusterPos i 
