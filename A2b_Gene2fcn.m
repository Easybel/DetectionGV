            % %                                  % %
            % % %                              % % %
            % % % % gene2fcn.m - by MonaIsa % % % %
            % % %                              % % %
            % %                                  % % 

% This script is for identifying gene functions of gene lists or
% MultiHitStat.mat - Matrices. The output cell BSU contains all genes
% {:,1} with their categories {:,2}; The output gcCat_tree contains the
% categories in column 1, the number of total hit (all genes) in column 2,
% the number of hits from your input gene list in column 3 and the
% subcategories in column 4

clear all; close all;
GC3PLOT=[5, 10, 30]; %[10, 12]; % If you like to have GC3 plots please type here ! e.g. [10, 12]  % To decide the number look at Cat2 (Variable)
%GC3PLOT = [1 : 33];% [16, 5, 22, 10, 9, 24, 23]; %[5, 10, 32]; % Vscy20 % [1, 3, 10, 13, 16, 22, 24, 32]; % Vnsc20                

savepath = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Plots_Default\'
savepath = 'J:\nullmodel\v2-noident\PLOTS\';

% Define your input genes - GeneOutput.mat:
%genes = 'J:\Wns_onBs166NCe\GeneOutput_WnsCy20_20200608_MF.mat';
%genes = 'J:\Ws_onBs166NCe\GeneOutput_Ws_Cy20_20200608.mat'; 
% genes = 'J:\sciebo\ResultsShared\Data_Summary\GeneOutput_Vns_cy20_woSample15.mat'; % Vnsc20
genes = 'J:\sciebo\ResultsShared\Data_Summary\GeneOutput_Vs_cy20_woSample9.mat'; % Vsc20
% Testing with simulated genes
%genes = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\NULLMODEL\Version 1 - No Ident only Adist average cluster number\20200709_TEST with simulated 14 replicates_Vnsc20\RUN2\SimulatedReplicates_GeneOutput.mat'; 

%
% % GeneOutput of simulated replicates (NullModel)
%
simgenes = 'J:\nullmodel\v2-noident\SimGeneOutput_Vnscy20_24july2020.mat';
simgenes = 'J:\nullmodel\v2-noident\SimGeneOutput_Vscy20_29july2020.mat';

% SubtiWiki gene categories in .txt-file
genecat = 'C:\Users\Mona\Documents\PhD_Mona Foerster\Bioinformatics\Dictionaries\genecategories\GeneCategories_2020_srt.txt';

% Would you like to have a Bonferroni correction for your p values?? % % % NEEDS TO BE REMOVED !
BONFERRONI = false; % true: with correction (default), false: without correction

% Saveme Abfrage
yes = 1; y = 1;
no = 0; n = 0;
saveme = input('Would you like to save the plots to savepath? ');
if saveme == 1
    run = input('Which experiment do you plot(e.g. Ws/Wns/Vs/Vnscy20)? ', 's');
end


         % %                                           % %
        % % %  Reading in SubtiWiki gene categories   % % %
         % %                                           % % 

fid = fopen(genecat); 
    inp = textscan(fid,'%s %s %s %s %s %s %s %s ', 'delimiter', '\t', 'HeaderLines', 1);
fclose(fid);
    gcBSU = inp{2}; %    gcGN = inp{3};
    gcCat1 = inp{4};
    gcCat2 = inp{5}; gcCat3 = inp{6};
clear inp

Cat1 = unique(gcCat1);
Cat2 = unique(gcCat2);
% And to which Cat1 do they belong ?
for i = 1 : length(Cat2)
    tmpmask = strcmp(Cat2{i,1}, gcCat2);
    parent = gcCat1(tmpmask);
    parent = unique(parent);
    Cat2{i,3} = parent;
end
Cat3 = unique(gcCat3);
% And to which Cat2 do they belong ?
for i = 1 : length(Cat3)
    tmpmask = strcmp(Cat3{i,1}, gcCat3);
    parent = gcCat2(tmpmask);
    parent = unique(parent);
    Cat3{i,3} = parent;
end

         % %                                             % %
        % % %  Reading in your genes of interest (GOI)  % % %
         % %                                             % % 

if strcmp(genes(end-2:end), 'txt') % THIS IS NOT USED ANYMORE
    % Loading BSU from txt - input file (BSU in 1. column!)
    ImpOpt = detectImportOptions(genes_ref);
    inp = readtable(genes, ImpOpt);
    BSU = table2cell(inp(:,1));
        
elseif strcmp(genes(end-2:end), 'mat') % THIS IS USED !
    % Loading BSU from geneoutput.mat
    inp = load(genes);
    inp = inp.GeneOutput;
    inpsamples = {inp(:).Sample}';
    inpsamples = unique(inpsamples);
%    inpsamples(10) = []; %% !!!!!!!!!!!!!!!!!!!THIS IS FOR EXCLUDING REPLICATE 10
    
end

             % %                                    % %
            % % %  Reading in the simulated genes  % % %
             % %                                    % % 
 
siminp = load(simgenes);
siminp = siminp.GeneOutput;
simsamples = {siminp(:).Sample};
simsamples = unique(simsamples);

         % %                                               % %
        % % %   Count categories of the simulated genes   % % %
         % %                                               % % 
     
    Cat1{1,2} = [];
    Cat2{1,2} = [];
    Cat3{1,2} = [];
for m = 1 : length(simsamples)
    cat1collect = []; cat2collect = [];cat3collect = [];
    tmp = simsamples{m};
    samplemask = strcmp(tmp, {siminp.Sample});
    BSU = [siminp(samplemask).BSU];
    
    for k = 1 : length(BSU)
        hitmask = strcmp(BSU(k), gcBSU);
        cat1collect = [cat1collect; gcCat1(hitmask)];
        cat2collect = [cat2collect; gcCat2(hitmask)];
        cat3collect = [cat3collect; gcCat3(hitmask)];
    end
    
    for i = 1 : size(Cat1,1)
        Cat1{i,2} = [Cat1{i,2} length(nonzeros(strcmp(Cat1(i), cat1collect)))];
    end
    for i = 1 : size(Cat2,1)
        Cat2{i,2} = [Cat2{i,2} length(nonzeros(strcmp(Cat2(i), cat2collect)))];
    end
        for i = 1 : size(Cat3,1)
        Cat3{i,2} = [Cat3{i,2} length(nonzeros(strcmp(Cat3(i), cat3collect)))];
    end
end
clear cat*collect

             % %                                     % %
            % % %  Plot the simulated distribution  % % %
             % %                                     % % 

figcounter = 1;

% % Plot Category 1 (Parent category)
figure(1); hold on;
set(gcf, 'renderer', 'painters', 'position', [10 10 1100 800]);

boxplot([Cat1{1,2}', Cat1{2,2}', Cat1{3,2}', Cat1{4,2}', Cat1{5,2}', Cat1{6,2}'], 'Labels', {Cat1{:,1}}, 'Notch', 'on', 'Color', 'k', 'Symbol', 'k+');
% Mean 
% plot([1:6], mean([Cat1{1,2}', Cat1{2,2}', Cat1{3,2}', Cat1{4,2}', Cat1{5,2}', Cat1{6,2}']), 'ko', 'LineWidth', 1.1, 'MarkerSize', 6);

set(gca,'xticklabelrotation',15);
set(gca, 'FontSize', 14)
 
% % Plot subset distributions - Category 2
for j = 1 : size(Cat1,1)
    figcounter = figcounter + 1;
    subset = j; % Which subset of Cat1 you want to plot ? 1: Cellular processes, 2: Groups of genes, ... (see Cat1{:,1})

        parentidx = strcmp(Cat1{subset,1},[Cat2{:,3}]);
        Cat2PlotSim{j} = vertcat(Cat2{parentidx,2})';
        %subsubidx = strcmp(Cat2{subsubset,1},[Cat3{:,3}]);
        Cat2PlotSimLabels{j} = {Cat2{parentidx, 1}};

    figure(figcounter); hold on;
    if j ~= 5 && j ~= 5
        set(gcf, 'renderer', 'painters', 'position', [50 50 max(size(Cat2PlotSim{j},2)*125, 950) 700]);
    elseif j == 5 || j == 3
        set(gcf, 'renderer', 'painters', 'position', [50 50 1100 700]);
    end
    boxplot(Cat2PlotSim{j}, 'Labels', Cat2PlotSimLabels{j}, 'Color', 'k', 'Notch', 'on', 'Symbol', 'k+');
    % Mean 
    % plot([1:size(Cat2PlotSim{j},2)], mean(Cat2PlotSim{j}), 'ko', 'LineWidth', 1.1, 'MarkerSize', 6);

    set(gca,'xticklabelrotation',15);
    set(gca, 'FontSize', 14)
    title(Cat1{subset,1});

    Cat2Plotmax(j) = max(Cat2PlotSim{j}(:));
end

% Clear two entries from the Cat3
    Cat3(1,:) = []; % These categories do not have a lower cat than Cat2
    Cat3(end-1,:) = []; % {'phosphorelay'} <-- belongs to 2 Cat2s: {'Regulation of gene expression';'Sporulation'}
% 
    Cat3{43,3} = 'Exponential and early post-exponential lifestyles'; % And this cat3 ('Genetic competence') belongs to 2 Cat2s ('Exp Lifestyles'+'Genetics') and we decide for only one 

    
% Plot category 3

for subsubset = GC3PLOT

    %subset = j; % Which subset of Cat1 you want to plot ? 1: Cellular processes, 2: Groups of genes, ... (see Cat1{:,1})

    parentidx = strcmp(Cat2{subsubset,1},[Cat3{:,3}]);
    Cat3PlotSim{subsubset} = vertcat(Cat3{parentidx,2})';
    Cat3PlotSimLabels{subsubset} = {Cat3{parentidx, 1}};
    if isempty(Cat3PlotSim{subsubset})
        continue
    end
    figcounter = figcounter + 1;
    figure(figcounter); hold on;
    set(gcf, 'renderer', 'painters', 'position', [50 50 max(size(Cat3PlotSim{subsubset},2)*125, 900) 700]);

    boxplot(Cat3PlotSim{subsubset}, 'Labels', Cat3PlotSimLabels{subsubset}, 'Color', 'k', 'Notch', 'on', 'Symbol', 'k+');
    % Mean 
    % plot([1:size(Cat3PlotSim{subsubset},2)], mean(Cat3PlotSim{subsubset}), 'ko', 'LineWidth', 1.1, 'MarkerSize', 6);

    set(gca,'xticklabelrotation',15);
    set(gca, 'FontSize', 14)
    title(Cat2{subsubset,1});

    Cat3Plotmax(subsubset) = max(Cat3PlotSim{subsubset}(:));
end

         % %                                   % %
        % % %  Count categories of your GOI   % % %
         % %                                   % % 
 
Cat1{1,3} = [];
Cat2{1,4} = [];
Cat3{1,4} = [];
for m = 1 : length(inpsamples)
    cat1collect = []; cat2collect = [];cat3collect = [];
    tmp = inpsamples{m};
    samplemask = strcmp(tmp, {inp.Sample});
    BSU = [inp(samplemask).BSU];
    
    for k = 1 : length(BSU)
        hitmask = strcmp(BSU(k), gcBSU);
        cat1collect = [cat1collect; gcCat1(hitmask)];
        cat2collect = [cat2collect; gcCat2(hitmask)];
        cat3collect = [cat3collect; gcCat3(hitmask)];
    end
    
    for i = 1 : size(Cat1,1)
        Cat1{i,3} = [Cat1{i,3} length(nonzeros(strcmp(Cat1(i), cat1collect)))];
    end
    for i = 1 : size(Cat2,1)
        Cat2{i,4} = [Cat2{i,4} length(nonzeros(strcmp(Cat2(i), cat2collect)))];
    end
    for i = 1 : size(Cat3,1)
        Cat3{i,4} = [Cat3{i,4} length(nonzeros(strcmp(Cat3(i), cat3collect)))];
    end
    
end

         % %                                      % %
        % % %  Plot the GOI hits and do ks-test  % % %
         % %                                      % % 
         
figcounter = 1;
figure(1); hold on;
Cat1Plot = vertcat(Cat1{:,3})';
Cat1PlotLabels = {Cat1{:,1}};
for i = 1 : size(Cat1,1)
    xvalue = ones(length(Cat1{i,3}),1)*i;
    plot(xvalue, Cat1{i,3}, 'bx', 'LineWidth', 1.5, 'MarkerSize', 6)
    plot(i, mean(Cat1{i,3}), '*', 'LineWidth', 2, 'MarkerSize', 12, 'Color', [0.6 0.6 1])
    plot(i, median(Cat1{i,3}), 'o', 'LineWidth', 2, 'MarkerSize', 12, 'Color', [0.6 0.6 1])
    
    [sigGC1, pGC1, ks2statGC1] = kstest2(Cat1{i,2}, Cat1{i,3});
    
    if BONFERRONI == true
    % Bonferroni - correction
    pGC1 = pGC1 * length(inpsamples);
    end
    
    if pGC1 < 0.05
        text(i-0.3, ceil(max(Cat1Plot(:))/100)*100 + 10, ['p = ', num2str(pGC1, '%.4f')], 'FontSize', 12, 'Color', 'r');
        disp(['GC1: ', Cat1{i,1}, ' is significantly deviating from your simulation data.']);
    elseif pGC1 >= 0.05
        text(i-0.3, ceil(max(Cat1Plot(:))/100)*100 + 10, ['p = ', num2str(pGC1, '%.4f')], 'FontSize', 12, 'Color', 'k');
    end
end

ylim([0 ceil(max(Cat1Plot(:))/100)*100 + 50]);
ylabel('Hits per replicate');

%boxplot(Cat1Plot, 'Labels', Cat1PlotLabels )

for j = 1 : size(Cat1,1)
    figcounter = figcounter + 1;
    subset = j;

    figure(figcounter); hold on;

    for i = 1 : size(Cat2,1)
        parentidx = strcmp(Cat1{subset,1},[Cat2{:,3}]);
        Cat2Plot = vertcat(Cat2{parentidx,4})';
        Cat2PlotLabels = {Cat2{parentidx, 1}};
    end

    for i = 1 : size(Cat2Plot,2)
        xvalue = ones(size(Cat2Plot,1),1)*i;
        plot(xvalue, Cat2Plot(:,i), 'bx', 'LineWidth', 1.5, 'MarkerSize', 6)
        plot(i, mean(Cat2Plot(:,i)), '*', 'LineWidth', 2, 'MarkerSize', 12, 'Color', [0.6 0.6 1])

        [sigGC2, pGC2, ks2statGC2] = kstest2(Cat2Plot(:,i), Cat2PlotSim{j}(:,i));
    
        if BONFERRONI == true
        % Bonferroni - correction
        pGC2 = pGC2 * length(inpsamples);
        end
    
        if pGC2 < 0.05
            text(i - 0.3, max(max(Cat2Plot(:)), Cat2Plotmax(j)) + 5, ['p = ', num2str(pGC2, '%.4f')], 'FontSize', 12, 'Color', 'r');
            disp(['GC2: ', Cat2PlotLabels{i}, ' is significantly deviating from your simulation data.']);
        elseif pGC2 >= 0.05
            text(i - 0.3, max(max(Cat2Plot(:)), Cat2Plotmax(j)) + 5, ['p = ', num2str(pGC2, '%.4f')], 'FontSize', 12, 'Color', 'k');
        end

    end

    ylim([0 max(max(Cat2Plot(:)), Cat2Plotmax(j))+ 10]);
    ylabel('Hits per replicate');  
end


for subsubset = GC3PLOT

    if isempty(Cat3PlotSim{subsubset})
        disp(['', Cat2{subsubset}, ' does not have any subcategories.']);
        continue
    end
    
    figcounter = figcounter + 1;
    figure(figcounter); hold on;

    for i = 1 : size(Cat3,1)
        parentidx = strcmp(Cat2{subsubset,1},[Cat3{:,3}]);
        Cat3Plot = vertcat(Cat3{parentidx,4})';
        Cat3PlotLabels = {Cat3{parentidx, 1}};
    end

    for i = 1 : size(Cat3Plot,2)

        xvalue = ones(size(Cat3Plot,1),1)*i;
        plot(xvalue, Cat3Plot(:,i), 'bx', 'LineWidth', 1.5, 'MarkerSize', 6)
        plot(i, mean(Cat3Plot(:,i)), '*', 'LineWidth', 2, 'MarkerSize', 12, 'Color', [0.6 0.6 1])

        [sigGC3, pGC3, ks2statGC3] = kstest2(Cat3Plot(:,i), Cat3PlotSim{subsubset}(:,i));
            
        if BONFERRONI == true
        % Bonferroni - correction
        pGC3 = pGC3 * length(inpsamples);
        end
    
        if pGC3 < 0.05
            text(i - 0.3, max(max(Cat3Plot(:)), Cat3Plotmax(subsubset)) + 5, ['p = ', num2str(pGC3, '%.4f')], 'FontSize', 12, 'Color', 'r');
            disp(['GC3: ', Cat3PlotLabels{i}, ' is significantly deviating from your simulation data.']);
        elseif pGC3 >= 0.05
            text(i - 0.3, max(max(Cat3Plot(:)), Cat3Plotmax(subsubset)) + 5, ['p = ', num2str(pGC3, '%.4f')], 'FontSize', 12, 'Color', 'k');
        end

    end

    ylim([0 max(max(Cat3Plot(:)), Cat3Plotmax(subsubset)) + 10]);
    ylabel('Hits per replicate');  
end
%% Saving the plots

if isempty(saveme) || saveme == 0 
    disp('Ok, your plots are not saved!');
elseif saveme == 1
    figure(1)
    print(gcf,  '-painters', '-dpng', [savepath, datestr(now, 'yyyymmdd'),'_', run, '_GC1']);
    print(gcf,  '-painters', '-depsc2', [savepath, datestr(now, 'yyyymmdd'),'_', run, '_GC1']);
    for j = 2 : figcounter
        figure(j);
        print(gcf,  '-painters', '-dpng', [savepath, datestr(now, 'yyyymmdd'),'_', run, '_', num2str(j)]);
        print(gcf,  '-painters', '-depsc2', [savepath, datestr(now, 'yyyymmdd'),'_', run, '_', num2str(j)]);
    end
else 
    disp('Ok, your plots are not saved!');
end
