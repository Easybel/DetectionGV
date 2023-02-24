%
%
% MULTIMAPPER.m by Mona
%
%
% To identify multimapper regions, you need to blastn the recipient on 
% itself and load the hittable(text) here. 
% After testing different filters and comparing with IGV (multimapper reads
% in white), we found the best values to be ident_min = 99.5 % and minlen
% (LIST2TABLES.m) = 300 to reproduce the bwa mem multimappers. 
%
% If you only like to include the donor on recip - multimappers as well, 
% please set 'include_donor' on 'y' (yes) (default: 'n' = no) and upload
% the blastn hitlist for donor on recip as well.
%
% Use LIST2TABLE.m to convert the output list to a table.

clear all; close all;

% Set paths and file names

blastn.path = 'H:\0_PhD\Bioinformatics\Multimapper_BLASTN\Bs166_NCe\';
blastn.recipname = 'Bs166NCe';
blastn.recip = 'HitTable_SbjBs166NCe_QryBs166NCe';
blastn.donorname = 'Bxx';
blastn.donor = 'HitTable_SbjBs166NCe_QryBsubW23';

output.path = blastn.path;
% output.path = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Multimapper\Batro\';

% Would you like to include the 'donor on recip' - multimappers(y/n)? (default:
% n)

include_donor = 'n';

% Parameters
ident_min = 99.5; % Set ident_min to the minimum identity (in percent) you want to allow 
sbj.size = 4227341; % Bs166_ass
%sbj.size = 4215607; % Bs166_NCe
%sbj.size = 3983858; % Bval
%sbj.size = 4168266; % Batro
%sbj.size = 3873116; % Geobacillus
%sbj.size = 4027680; % BsubW23

% how long should a multimapper at least be?
minLen = 40;


%%
% Read in data
% 

% BlastN Hittable of blasting recipient on recipient
fid = fopen([blastn.path, blastn.recip, '.txt']);
    input1 = textscan(fid, '%s %s %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 8);
    sbj.start1 = input1{9};
    sbj.end1 = input1{10};
    ident1 = input1{3};
fclose(fid);

% BlastN Hittable of blasting donor on recipient
if include_donor == 'y'
    fid = fopen([blastn.path, blastn.donor, '.txt']);
        input2 = textscan(fid, '%s %s %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 8);
        sbj.start2 = input2{9};
        sbj.end2 = input2{10};
        ident2 = input2{3};
    fclose(fid);
end


% Filter out parts where the identity is too low, if you like to do so, if
% not please set ident_min = 0 at the start of the script!
lowident_filt1 = find(ident1>=ident_min);
sbj.start1 = sbj.start1(lowident_filt1);
sbj.end1 = sbj.end1(lowident_filt1);

if include_donor == 'y'
    lowident_filt2 = find(ident2>=ident_min);
    sbj.start2 = sbj.start2(lowident_filt2);
    sbj.end2 = sbj.end2(lowident_filt2);
end

clear ident_min
clear lowident_filt1 lowident_filt2;

% How many queries cover which position on the recipient
% 1. for recip on recip:
hit.pos1 = zeros(sbj.size,1);   % hit.pos1 contains the number of hits per position of the recipient on itself
for i=1:length(sbj.start1)
    if sbj.start1(i) >= sbj.end1(i)
        hit.pos1(sbj.end1(i):sbj.start1(i)) = hit.pos1(sbj.end1(i):sbj.start1(i))+1;
    end
    if sbj.start1(i) < sbj.end1(i)
        hit.pos1(sbj.start1(i):sbj.end1(i)) = hit.pos1(sbj.start1(i):sbj.end1(i))+1;
    end
end

% 2. for donor on recip:
if include_donor == 'y'
    hit.pos2 = zeros(sbj.size,1);   % hit.pos1 contains the number of hits per position of the donor on the recip
    for i=1:length(sbj.start2)
        if sbj.start2(i) >= sbj.end2(i)
            hit.pos2(sbj.end2(i):sbj.start2(i)) = hit.pos2(sbj.end2(i):sbj.start2(i))+1;
        end
        if sbj.start1(i) < sbj.end1(i)
            hit.pos2(sbj.start2(i):sbj.end2(i)) = hit.pos2(sbj.start2(i):sbj.end2(i))+1;
        end
    end
end
%% Identifying multimapper regions

multimapper1 = find(hit.pos1>1);
if include_donor == 'y'
    multimapper2 = find(hit.pos2>1);
end

if include_donor == 'y'
    multimapper = vertcat(multimapper1, multimapper2);
    disp('Your list includes "recip on recip" - and "donor on recip" - multimappers!');
elseif include_donor == 'n'
    multimapper = multimapper1;
    disp('Your list only includes "recip on recip" - multimappers!');
else 
    error('Please, set exclude_donor either on y(es) or n(o)!');
end

% Sort out doubled positions
multimapper = unique(multimapper);

% How much of the genome is multimapper region ? 
perc_multim = length(multimapper) / sbj.size * 100;
disp([num2str(perc_multim, '%0.3f'), ' % of the genome are multimapper regions.']);
disp(['That means ', num2str(length(multimapper), '%i'), ' basepairs are affected.']);

%% Write to file (Use LIST2TABLE.m to make a table out of this list)
if include_donor == 'y'
    output.name = [output.path, blastn.recipname, '_', blastn.donorname, '_donorincl_mm_list.txt'];
elseif include_donor == 'n'
    output.name = [output.path, blastn.recipname, '_mm_list.txt'];
end 

outfasta = fopen(output.name, 'w');
fprintf(outfasta, '%d\n', multimapper);
fclose(outfasta); 


%% convert the list of positions to a list with start and end

ListIdx_start = find([multimapper - [1;multimapper(1:end-1)]] ~=1);
ListIdx_ende = [ListIdx_start(2:end) - 1; numel(multimapper)];

ListRegion(:,1) = multimapper(ListIdx_start);
ListRegion(:,2) = multimapper(ListIdx_ende);
ListLength =  ListRegion(:,2) - ListRegion(:,1) +1;

Listmask = ListLength >= minLen;

output.nameList = [output.path, blastn.recipname, '_mm.txt'];
outList = fopen(output.nameList, 'w');
fprintf(outList, '%d %d %d\n', [ListRegion(Listmask,:) ListLength(Listmask)]');
fclose(outList); 
