%
% JJPsBoxplot.m
%
clear all; close all;

%Here you can add as many Cluster information as you want
CNPname{1} = 'J:\Vns_onBs166NCe\CNPSummary.mat'; 
%CNPname{1} = 'J:\Ws_onBs166NCe\CNPSummary_Ws_Cy20_23042020_MeY.mat';
%CNPname{2} = 'J:\Ws_onBs166NCe\CNPSummary_Ws_Cy10_21042020_MeY.mat';

% Which samples you want to plot? - You can use one or more data sets
%dataset{1} = {'Ws0610','Ws0710','Ws0810','Ws0910', 'Ws1010', 'Ws1210'}; % data set 1
%dataset{2} = {'Ws0620','Ws0720','Ws0820','Ws0920', 'Ws1020', 'Ws1220'}; % data set 2
dataset{1} = {'Vns0110', 'Vns0210','Vns0310', 'Vns0410','Vns0510', 'Vns0610','Vns0710','Vns0810', 'Vns0910', 'Vns1010', 'Vns1110', 'Vns1210', 'Vns1310','Vns1410','Vns1510'}; % data set 3

sampleP = [];
for i = 1 : length(dataset)
    sampleP = [sampleP dataset{i}];
end

%% Calculations ...

% Read in data
Cluster = [];
for i = 1 : length(CNPname)
    Data = load(CNPname{i});
    Cluster = [Cluster; Data.CNPSummary];
end
clear Data;

% Are the samples you want to plot in the CNPSummaries you loaded? - log
% And, where in CNPname are the samples you want to plot? - pos 
pos = [];
for i = 1 : length(dataset)
[log, pos_tmp] = ismember(dataset{i}, {Cluster(:).Samples});
ERR1 = find(log == 0);

% Make a list of the sample positions in the Cluster struct (in the right
% order)
pos = [pos pos_tmp];

if ~isempty(ERR1)
    error(['At least one sample of your data set ', num2str(i), ' is not in the loaded CNPSummaries.']);
end
end
clear ERR1 pos_tmp log

% Collect only the clusters for this plot in CNP_plot (and put them in 
% the right order)
CNP_plot = Cluster(pos);

% For better use, split the CNPsummaries:
% Create an Adist cell
Adist = {CNP_plot(:).Adist};

% % Boxplots need special arrays:
% Create an array that can be used for boxplots (all values in one column)
adist_all = [];
% ...and an array with the name and number of the sample the values belongs 
% to 
adist_repl = [];
adist_repl_int = [];

for i = 1 : length(CNP_plot)
    % Collect all adist values (of all data sets)
    adist_all = [adist_all; Adist{i}(:,1)];
    if length(sampleP{i}) == 6
        tmp_name = [char(sampleP(i)) '  '];
    elseif length(sampleP{i}) == 7
        tmp_name = [char(sampleP(i)) ' '];
    elseif length(sampleP{i}) == 8
        tmp_name = [char(sampleP(i))];
    end
    % .. and their names
    for j = 1 : length(Adist{i}(:,1))
        adist_repl = [adist_repl; tmp_name];
    end
    % .. and the sample number (1 to length(sampleP))
    adist_repl_int = [adist_repl_int; ones(length(Adist{i}),1)*i];
end
% Now, the arrays for the boxplot are done.


%% Plotting ...
%
% Figure 1 : Box-whisker plot of mean integration length for each sample +
%            Percentage of genome replaced (JJP)
% 

p1 = figure(1); hold on;
set(p1, 'renderer', 'painters', 'position', [20 20 115*length(sampleP) 500],'paperpositionmode','auto'); 
set(gca, 'FontSize', 16);

% Boxplot
boxplot(adist_all, adist_repl, 'whisker', 1, 'plotstyle', 'traditional');

% Plotting every Adist point per sample, to not have all at the same x
% position, shift them with the rand function
plot(adist_repl_int + 0.4 .* (rand(length(adist_repl_int),1) - 0.5), adist_all, 'ko');
ylabel('Mean integrated segment length')

% On the right site, plot the percentage of genome replaced in green
yyaxis right
set(gca, 'YColor', [0 0.79 0.4]); 
plot([CNP_plot(:).Transfer], '*', 'color', [0 0.79 0.4], 'LineWidth', 2, 'MarkerSize', 10);
ylabel('Perc. genome replaced (green *)');
ylim([0 round(max([CNP_plot(:).Transfer])+2)]);
