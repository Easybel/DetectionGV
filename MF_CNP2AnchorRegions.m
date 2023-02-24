 % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
 % % % % % % % % BE CAREFUL! THIS IS JUST A BETA VERSION % % % % % % % % 
 % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
 
   %                                                               %
   % %                                                           % %
   % % %   Identity of start and end of our detected clusters  % % %
   % %                                                           % %
   %                                                               %
clearvars; close all

%
% %
% % % Variables and paths
% %
%

% Define your length of interest (LOI) for the start and edge segment
LOI = 26; %[bp]

recipsize = 4215607;
savepath = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Plots_Default\';
mmlist = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Dictionaries\mmlist\Bs166NCe_mm.txt';

% % Tell us more about the samples you would like to test ... 
% % - CORRECTED wants to know if you use CNPSummaries with the latest Mdist correction (yes) or an old one (no) 

% CNPname = 'J:\Vns_onBs166NCe\20200727_Vns_Cy10_Cy20_CNPSummary_MF.mat'; 
% dataset = [21 : 34];
% setname = 'Vnsc20'
% donor = 'v';
% CORRECTED = 'yes'

% CNPname = 'J:\Vns_onBs166NCe\20200727_Vns_Cy10_Cy20_CNPSummary_MF.mat'; 
% dataset = [1 : 15];
% setname = 'Vnsc10'
% donor = 'v';
% CORRECTED = 'yes'

CNPname = 'J:\sciebo\ResultsShared\Data_Summary\CNPSummary_Vscy10_20_DP50.mat';
dataset = 21:35;
setname = 'Vsc20';
donor = 'v';
CORRECTED = 'no'

% CNPname = 'J:\Ws_onBs166NCe\CNPSummary_Ws_CY10_CY20_19052020_MeY.mat';
% dataset = 21:35;
% setname = 'Wsc20';
% donor = 'w';
% CORRECTED = 'no'

% CNPname = 'J:\sciebo\ResultsShared\Data_Summary\CNPSummary_Wns_CY10_CY20_14052020_MF.mat';
% dataset = 17:22;
% setname = 'Wnsc20';
% donor = 'w';
% CORRECTED = 'no'

if isempty(donor)
    donor = input('What donor is used in your experiment? w: BsubW23, v: Bval, a: Batro ', 's');
end
if donor == 'w'     % W23 donor
    masterlist = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Dictionaries\masterlists\BsubW23_2Bs166NCe_DP50_ml.txt';
    accgenome = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Dictionaries\accgens\BsubW23_2Bs166NCe_DP50_acc.txt';
elseif donor == 'v' % Bval donor
    masterlist = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Dictionaries\masterlists\Bval_2Bs166NCe_DP50_ml.txt';
    accgenome = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Dictionaries\accgens\Bval_2Bs166NCe_DP50_acc.txt';
elseif donor == 'a' % Batro donor
    error('You need to renew master list and accessory genome!');
else 
    error('You need to specify your donor for the following analysis! Please try again.')
end

%
% %
% % % Read in - Section ...
% %
%

Data = load(CNPname);
CNPSummary = Data.CNPSummary(dataset);
clear Data;

fid = fopen(masterlist);
imp = textscan(fid,'%f %s %s');
fclose(fid);
ml.pos = imp{1};
ml.atcg = imp{3};
clear imp

fid = fopen(accgenome);
imp = textscan(fid, '%f %f %f', 'delimiter', '\t');
acc.start = imp{1};
acc.edge = imp{2};
fclose(fid);
clear imp

fid = fopen(mmlist);
imp = textscan(fid, '%f %f %f', 'delimiter', '\t');
mm.start = imp{1};
mm.edge = imp{2};
fclose(fid);
clear imp

%
% % 
% % % Calculations ...
% %
%

% Calculate local identity 
% (1: identical allele, 0: alternating allele)

% Create array with ones for every position on the genome
id = ones(1, recipsize);

% All positions on the master list differ from recipient alleles
id(ml.pos) = 0;

% All positions within acc or mm(??) regions differ from recipient alleles
for i = 1 : length(acc.start)
    id([acc.start(i):acc.edge(i)]) = 0;
end
for i = 1 : length(mm.start)              %% Do we need to set mm regions to 0 identity ?
    id([mm.start(i):mm.edge(i)]) = 0;
end

% Use moving mean to calculate the average identity of a LOI long region
% at position pos : locid = movmean(id, [0 seglen - 1])
locid = movmean(id, [0 LOI - 1]);

% How many different start positions on the whole genome do we have with 
% 'LOI' identical alleles in a row?
MEPS_tot = length(locid(locid==1))

% Find all start and edge points of the experimental data
Mdist = vertcat(CNPSummary(:).Mdist);
% Filter out clusters with mdist < LOI
mask = (Mdist(:,1) >= LOI);
Mdist_start = Mdist(mask,2);
Mdist_edge = Mdist(mask,3);
clear Mdist;

corrtest = 0;
while corrtest == 0
if strcmp(CORRECTED, 'yes')
    % for the corrected CNPSummaries:
    id_start = locid(Mdist_start);
    id_edge = locid(Mdist_edge - LOI + 1);
    corrtest = 1;
elseif strcmp(CORRECTED, 'no')
    % for the uncorrected CNPSummaries:
    id_start = locid(Mdist_start + 1);
    id_edge = locid(Mdist_edge - LOI); % + 1);
    corrtest = 1;
else 
    CORRECTED = input('Is your CNPSummary correct?? - yes/no', 's');
end
end

% Decide for the higher identity - Potential anchor region
id_max = max(id_start, id_edge);
num_max = id_max * LOI;
id_min = min(id_start, id_edge);
num_min = id_min * LOI;

%
% %
% % % Plotting ....
% %
%

% Set histogram parameters depending on LOI value
if LOI <= 70
    BINWIDTH = 1; BW2val = 1/LOI; BW2 = BW2val;
elseif LOI > 70 && LOI <= 140
    BINWIDTH = 2; BW2val = 1/LOI; BW2 = BW2val * 2;
elseif LOI > 140
    BINWIDTH = 5; BW2val = 1/LOI; BW2 = BW2val * 4;
end

figure(1); hold on;
title([setname, ' - Distribution of the front and edge segments']);
set(gcf, 'renderer', 'painters')
set(gca, 'FontSize', 12)
histogram(id_start*LOI, 'BinEdges', [BINWIDTH/2: BINWIDTH: LOI+BINWIDTH/2], 'FaceColor', [0.6 0.6 0.6],'EdgeAlpha', 0)
histogram(id_edge*LOI, 'BinEdges', [BINWIDTH/2: BINWIDTH: LOI+BINWIDTH/2], 'FaceAlpha', 0, 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.4)
xlim([min([id_start id_edge]*LOI)-0.5 LOI+BINWIDTH/2])
legend(['Mdist Front\newlineMean = ', num2str(mean(id_start*LOI), '%0.1f'), ' \pm ', num2str(std(id_start*LOI), '%0.1f')], ['Mdist Back\newlineMean = ', num2str(mean(id_edge*LOI), '%0.1f'), ' \pm ', num2str(std(id_edge*LOI), '%0.1f') ], 'Location', 'nw')
xlabel(['Number of fitting alleles in the first/last ', num2str(LOI), ' bps'])
ylabel('Probability')

figure(24); hold on;
title([setname, ' - Total number of fitting alleles']);
set(gcf, 'renderer', 'painters')
set(gca, 'FontSize', 12)
h1 = histogram(num_max,'BinEdges', [BINWIDTH/2: BINWIDTH: LOI+BINWIDTH/2], 'Normalization', 'probability', 'FaceAlpha', 0 ,'EdgeColor', [0 0 0], 'EdgeAlpha', 0.9, 'LineWidth', 1.6)
h2 = histogram(num_min,'BinEdges', [BINWIDTH/2: BINWIDTH: LOI+BINWIDTH/2], 'Normalization', 'probability', 'FaceColor', [0 0 0], 'FaceAlpha', 0.5, 'LineWidth', 1.6, 'EdgeAlpha', 0)
ylabel('Probability')
xlabel(['Number of fitting alleles in the first/last ', num2str(LOI), ' bps'])
legend('Potential anchor region', 'Potential end of recombination', 'Location','nw')
xlim([min(num_min) LOI+BINWIDTH/2])

figure(12); hold on;
title([setname, ' - All start and edge segments together']);
set(gcf, 'renderer', 'painters')
set(gca, 'FontSize', 12)
h3a = histogram(locid(locid>0), 'BinEdges', [BW2val/2: BW2: 1+BW2val/2], 'Normalization', 'probability', 'FaceColor', [0 0.6 0.3], 'FaceAlpha', 0.4, 'EdgeAlpha', 0)
h6 = histogram([id_min id_max], 'BinEdges', [BW2val/2: BW2: 1+BW2val/2], 'Normalization', 'probability', 'FaceAlpha', 0, 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1.6)
ylabel('Probability')
xlim([min([id_start id_edge])-0.1 1+BW2val/2])
xlabel(['Sequence identity of ', num2str(LOI), ' bp long segments'])
legend(['Average local ident (>0)\newlineMean = ', num2str(mean(locid(locid>0)), '%0.2f'), ' \pm ', num2str(std(locid(locid>0)),'%0.2f')], ['Recombination start + end \newlineMean = ', num2str(mean([id_min id_max]), '%0.2f'), ' \pm ', num2str(std([id_min id_max]), '%0.2f')], 'Location', 'nw')
% % If you like to test the distribution of the possible sequence identity values
% htest = histogram([1/LOI: 1/LOI: 1],'BinEdges',  [BW2cal/2: BW2: 1+BW2cal/2], 'Normalization', 'probability', 'FaceColor', [0 0 0]);

figure(22); hold on;
title(setname);
set(gcf, 'renderer', 'painters')
set(gca, 'FontSize', 12)
h3 = histogram(locid(locid>0), 'BinEdges', [BW2val/2: BW2: 1+BW2val/2], 'Normalization', 'probability', 'FaceColor', [0 0.6 0.3], 'FaceAlpha', 0.4, 'EdgeAlpha', 0)
h4 = histogram(id_min, 'BinEdges', [BW2val/2: BW2: 1+BW2val/2], 'Normalization', 'probability', 'FaceAlpha', 0, 'EdgeColor', [0.6 0.6 0.6], 'LineWidth', 1.6)
h5 = histogram(id_max, 'BinEdges', [BW2val/2: BW2: 1+BW2val/2], 'Normalization', 'probability', 'FaceAlpha', 0, 'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 1.6)
ylabel('Probability')
xlim([min([id_start id_edge])-0.1 1+BW2val/2])
xlabel(['Sequence identity of the first/last ', num2str(LOI), ' bps'])
legend(['Average local ident (>0)\newlineMean = ', num2str(mean(locid(locid>0))*100, '%0.1f'), ' \pm ', num2str(std(locid(locid>0))*100,'%0.1f'), ' %'], ['Potential end of recombination \newlineMean = ', num2str(mean(id_min)*100, '%0.1f'), ' \pm ', num2str(std(id_min)*100, '%0.1f'), ' %'], ['Potential anchor region \newlineMean = ', num2str(mean(id_max)*100, '%0.1f'), ' \pm ', num2str(std(id_max)*100, '%0.1f'), ' %'], 'Location', 'nw')

%
% %
% % % Saving the figures ...
% %
%
figs = [1 12 22 24];

printOpt = 0; % printOpt = 1;
if printOpt == 1;
    for i = 1 : length(figs)
        figure(figs(i));
        print(gcf,  '-painters', '-dpng', [savepath, datestr(now, 'yyyymmdd'),'_IDstart_', setname, '_LOI', num2str(LOI), '_', num2str(i)]);
    end
end

%
% %
% % % MEPS - Counting the consecutively fitting alleles no both ends
% %
%

% Create a list of all zero positions (mismatches) in the id array
mismatches = find(id == 0); % contains mlSNPs, accgen + mm regions

% Count MEPS_start and MEPS_edge  % % Another way would be to substract cdist and mdist ... = MEPS !!
for i = 1 : length(Mdist_start)
    % Find the smallest mismatch greater than mdist_start and calculate distance:
    MEPS_start(i) = min(mismatches(mismatches > Mdist_start(i))) - Mdist_start(i);  
    % Find the closest mismatch smaller than mdist_edge and calculate distance:
    MEPS_edge(i) = Mdist_edge(i) - max(mismatches(mismatches < Mdist_edge(i)));
end

MEPS = max(MEPS_start, MEPS_edge); % Decide for the greater value as potential start of recombination

% What do we learn about the MEPS of this experiment?
MEPS_min = min(MEPS);
MEPS_10 = length(nonzeros(MEPS<=10));
MEPS_15 = length(nonzeros(MEPS<=15));
MEPS_20 = length(nonzeros(MEPS<=20));
disp(['Smallest MEPS: ', num2str(MEPS_min)]);
disp(['MEPS smaller 10: ', num2str(MEPS_10), ', this is ', num2str(MEPS_10/length(Mdist_start)*100), ' %']);
disp(['MEPS smaller 15: ', num2str(MEPS_15), ', this is ', num2str(MEPS_15/length(Mdist_start)*100),' %']);
disp(['MEPS smaller 20: ', num2str(MEPS_20), ', this is ', num2str(MEPS_20/length(Mdist_start)*100), ' %']);

%
% %
% % % Plotting ....
% %
%

figure(2); hold on; 
set(gca, 'FontSize', 12)
set(gcf, 'Renderer', 'Painters')
histogram(MEPS, 'BinEdges', [0.5:1:round(max(MEPS))+0.5], 'Normalization', 'Count', 'EdgeAlpha', 0, 'LineWidth', 1.4, 'FaceColor', [0.3 0.3 0.3])
xlim([0 30])
legend(['Total number\newlineof clusters = ', num2str(length(Mdist_start))], 'Location', 'Nw');

figure(3); hold on; 
set(gca, 'FontSize', 12)
set(gcf, 'Renderer', 'Painters')
histogram(MEPS, 'BinWidth', 15, 'Normalization', 'Prob', 'EdgeAlpha', 0, 'LineWidth', 1.4, 'FaceColor', [0.3 0.3 0.3])
xlim([0 120])
legend(['Total number\newlineof clusters = ', num2str(length(Mdist_start))], 'Location', 'ne');
xlabel('Number of consecutively fitting SNPs in the potential anchor region');
ylabel('Probability');

%
% %
% % % 'Difftest' - What is the distance between two mlSNPs? (To answer why we see the bias in fig 2!)
% %
%

difftest = diff(ml.pos);
figure(100); hold on;
set(gca, 'FontSize', 12);
set(gcf, 'Renderer', 'painters');

histogram(difftest, 'BinEdges', 0.5:1:max(difftest)+0.5)
xlabel('Number of homologous alleles between mlSNPs');
xlim([0.5 60.5])
legend('A multiple of 3 homologous SNPs\newlinebetween mlSNPs seems to be favored');
