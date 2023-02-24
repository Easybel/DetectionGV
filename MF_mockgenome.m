%
% Create your own MOCKGENOME - by Mona
%
% This script reads in the dictionary of your recipient and your donor species and  
% the list of changes (start/end recip, start/end donor, ident) you want to have in 
% your mock genome. The output is: 
%   1.  the new mock fasta (you can then cut into reads with wgsim) and
%   2.  a .mat-matrix containig start/end recip, start/end donor, length recip, ident 
%       of the integrated pieces (for further use in matlab)

% At the moment: Optimized for v003 !!!

close all; clear all

% Define paths and variables
% Where do I find the reference dictionary you want to change? ('subject')
inp.path1 = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Dictionaries\';
inp.name1 = 'BsubNC_000964wt';
inp.headerlinelength1 = 1; % How many lines are belonging to the header?

% Where do I find the dictionary the new nt sequences come from? ('query')
inp.path2 = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Dictionaries\';
%inp.name2 = 'W23wt';
inp.name2 = 'Bvallismortis_DSM11031';
inp.headerlinelength2 = 1; % How many lines are belonging to the header?

% Where do you want me to print the outputs and what is the name of the project?
output.path = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Mockgenomes\Version003\';
output.name = 'Bs166_mock_v003';

% Please create a list of changes with 4 columns: query start, query end,
% subject start, subject end, identity (You can copy/paste it from the
% BLASTN HitTable (text) or add more by hand
change.path = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Mockgenomes\Version003\';
change.list = 'changelist_v003.txt';
change.headerlinelength = 1; % How many lines are belonging to the header?

% Read in subject fasta and delete header and '?' 
% The subject fasta is the fasta where some new parts from the query fasta
% are integrated - it gives you the basis of the mockgenome
fid = fopen([inp.path1, inp.name1, '.fasta'], 'r');
for i=1:inp.headerlinelength1  % skipping head line
    fgetl(fid);
end
nt = fread(fid, '*char');       % reading in the nucleotides nt of the dictionary
fclose(fid);

deleteme = find(nt~=10);        % delete the "enter" ('?': ascii = 10) from file
nt = nt(deleteme);
inp.length1 = length(nt);      % Parameter to check if everything went right
clear i deleteme ans fid;

% Read in query fasta and delete header and '?'
% The query fasta is the fasta where we cut some parts out and then
% integrate them into the subject fasta
fid = fopen([inp.path2, inp.name2, '.fasta'], 'r');
for i=1:inp.headerlinelength2  % skipping head line
    fgetl(fid);
end
nt_qry = fread(fid, '*char');       % reading in the nucleotides nt of the dictionary
fclose(fid);

deleteme = find(nt_qry~=10);        % delete the "enter" ('?': ascii = 10) from file
nt_qry = nt_qry(deleteme);
inp.length2 = length(nt_qry);       % Parameter to check if everything went right
clear i deleteme ans fid;

% Here, the GC content of the recipient fasta is calculated: 
G = length(find(nt=='G'));
A = length(find(nt=='A'));
T = length(find(nt=='T'));
C = length(find(nt=='C'));
GC = (G + C)/length(nt)*100;
clear G A T C

% Read in your changelist  
fid = fopen([change.path, change.list]);
changelist = textscan(fid, '%f %f %f %f %f %f %f %f %f', 'HeaderLines', 1);      % reading in the change list (qry.start/end, sbj.start/end, ident)
fclose(fid);

change.sbj(:,1) = changelist{1,1};  % Start position of subject
change.sbj(:,2) = changelist{1,2};  % End position of subject
change.qry(:,1) = changelist{1,3};  % Start position of query
change.qry(:,2) = changelist{1,4};  % End position of query
change.ident = changelist{1,9};     % Identity between the 2 pieces
clear changelist

% In some cases, subject and query pieces will not run in the same
% direction. On the blastn hittable, the subject then starts with a higher
% number than it ends with. To be able to sort the changelist in this
% script, we first need to switch subject starts and ends where this is the
% case. To keep the information, the query starts and ends then have to be
% switched aswell:
for i=1:length(change.sbj)
    if change.sbj(i,1) > change.sbj(i,2)
        temp = change.sbj(i,1);
        change.sbj(i,1) = change.sbj(i,2);
        change.sbj(i,2) = temp;
        temp = change.qry(i,1);
        change.qry(i,1) = change.qry(i,2);
        change.qry(i,2) = temp;
    end 
end

% sorting the change list because we integrate the new sequences beginning with
% the last one goint to the first
[x,sort_idx] = sort(change.sbj(:,1), 'descend');

change_sort.sbj = change.sbj(sort_idx, :);
change_sort.qry = change.qry(sort_idx, :);
change_sort.ident = change.ident(sort_idx);

clear x sort_idx

%
% See if you can find overlaps of the subject sequences, if so: delete
%
% If the start of the next piece is lower than the end of the piece before: overlaps > 0
overlaps = change_sort.sbj(2:end,2) - change_sort.sbj(1:end-1,1);  
% Finds the indices of overlapping pieces; Add one to the index so that the 
% second piece of the two overlapping pieces is hit
overlaps_idx = find(overlaps > 0) + 1;
% Create a mask to filter out the second hits
overlap_mask = ~ismember(change_sort.sbj(:,1),change_sort.sbj(overlaps_idx,1)); 

change_sort.sbj = change_sort.sbj(overlap_mask, :);
change_sort.qry = change_sort.qry(overlap_mask,:);
change_sort.ident = change_sort.ident(overlap_mask);

% Calculate the lengths
change_sort.qry(:,3) = change_sort.qry(:,2) - change_sort.qry(:,1) + 1;    % calculating the length of the integrated piece
change_sort.sbj(:,3) = change_sort.sbj(:,2) - change_sort.sbj(:,1) + 1;    % calculating length on subject (neg. value: from + to -)

clear overlap_mask overlaps overlaps_idx

% Create a new dictionary including the changes from the changelist
neg_qry = find (change_sort.qry(:,3) < 0); % Here the length is only right for the positive values . The negative values are 2 too short (e.g. -220 has a length of 222)
neg_sbj = find (change_sort.sbj(:,3) < 0); % Should be empty !
if ~isempty(neg_sbj)
    error('Ohoh! Something went wrong here! We cannot start integrating the new pieces. Check the script, girl')
end 

% Create the complement of the qry.fasta
nt_qry_compl = char;
idxA = find(nt_qry == 'A');
nt_qry_compl(idxA,1) = 'T';
idxT = find(nt_qry == 'T');
nt_qry_compl(idxT,1) = 'A';
idxC = find(nt_qry == 'C');
nt_qry_compl(idxC,1) = 'G';
idxG = find(nt_qry == 'G');
nt_qry_compl(idxG,1) = 'C';
if length(idxA) + length(idxT) + length(idxC) + length(idxG) ~= length(nt_qry)
    error('The created complementary query dictionary is not complete! Maybe there are other letters than ATCG in the dict (N?)');
end


nt_old = length(nt)
for i=1:length(change_sort.sbj(:,1))
    if change_sort.qry(i,3)>0
        nt = [nt(1:change_sort.sbj(i,1)-1)' nt_qry(change_sort.qry(i,1):change_sort.qry(i,2))' nt(change_sort.sbj(i,2)+1:end)']';
    elseif change_sort.qry(i,3)<0
        % for plus/minus - replacements you need to integrate the reverse
        % complement !:
        nt = [nt(1:change_sort.sbj(i,1)-1)' flipud(nt_qry_compl(change_sort.qry(i,2):change_sort.qry(i,1)))' nt(change_sort.sbj(i,2)+1:end)']';
    else error('You want to integrate pieces with length 0. Is this really what you want?');
        
    end
end

nt_new = length(nt)


%% Create output in .fasta - format
yes = 1; no = 0;
saveme = input(['Would you like to save the new fasta (yes/no)? ']);
if isempty(saveme) || saveme == 1;
    outfasta = fopen([output.path, output.name, '.fasta'], 'w');
    fprintf(outfasta, ['>', inp.name1, ' with ', inp.name2, '-pieces_MOCKFASTA\n']);
    fprintf(outfasta, '%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c\n', nt);
    fclose(outfasta); 
    % a .mat - matrix called changes_table is saved, it contains the start and
    % end of subject, start and end of query, length of subject piece and the
    % identity
    changes_table = [change_sort.sbj(:,1) change_sort.sbj(:,2) change_sort.qry(:,1) change_sort.qry(:,2) change_sort.sbj(:,3) change_sort.ident];
    changes_table = flipud(changes_table);
    save([output.path, output.name, '_changes.mat'],'changes_table');
else
    disp('Ok, the mock genome is not saved.')
end

% Tidy up
clear outfasta fid ans i temp neg_sbj nt_qry_flipped 
clear nt_qry output inp 
