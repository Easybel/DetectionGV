%
%
% PLOTSCRIPT.m - by MonaIsa
%

% % This script is for all the small plots you like to have for your
% sequencing data.
% Input: 
% -- Cluster information: CNPSummary of your replicates

% Order of this script: 1) all settings (Here, it is your turn!)
%                       2) all calculations
%                       3) all figures
clear all; close all;

% Where would you like to save the plots ?
savepath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/DNASeq/plots/";
basePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/DNASeq/2a_Cluster/";  
% % Samples

% Load your CNPSummaries
CNPname{1} = basePath + "20221229_Wns_CNPSummary.mat"; 
CNPname{2} = basePath + "20221229_Vns_CNPSummary_wo08.mat"; 
CNPname{3} = basePath + "20221229_BAns_CNPSummary.mat"; 

% Which samples you want to plot? - You can use one or more data sets
% Wns
% dataset{1} = ["Wns" + string([11:12 14:19]) + "10"]
dataset{1} = ["Wns" + string([11:12 14:19]) + "20.5"]
% Vns
% dataset{1} = ["Vns0" + string([1:8]) + "10"]
dataset{2} = ["Vns0" + string([1 3:4]) + "20" "Vns0" + string([2 5:7]) + "20.5"]
% BAns
% dataset{1} = ["BAns0" + string([1:4 6:8]) + "10"]
dataset{3} = ["BAns0" + string([1:4 6:8]) + "20"]


% Would you like to name your datasets ?
% Wns
% ds_name{1} = 'Wns Cy10';
ds_name{1} = 'Bspiz Cy20';
% Vns
% ds_name{1} = 'Vns Cy10';
ds_name{2} = 'Bval Cy20';
% BAns
% ds_name{1} = 'BAns Cy10';
ds_name{3} = 'Batro Cy20';

% Create a cell with all the samples in your datasets (sampleP)
sampleP = [];
for i = 1 : length(dataset)
    sampleP = [sampleP dataset{i}];
end

% % Colors

% Create your own colormap
clrmp = [0.1 0.1 0.1; 0 0.4 1; 0 0.79 0.61; 0.9 0 0; 0 0 0.6; 0 0 0.8; 0 0 1; 0 0.2 1; 0 0.4 1; 0 0.6 1];
cmp = repmat([189 38 38; 0 76 255; 81 189 81; 232 163 23; 0 0 0;  235 235 84; 130 218 250; 247 134 247]/255, 70,1);

% Or choose an existing one:
% Color per dataset
clr_ds = clrmp; %%summer(length(dataset)); 
% Color for every replicate on its own
clr_samp = jet(length(sampleP));

% % Identity histogram
% Do you want to have a minimum segments length for the identity histogram?
min_len = 100; % [bp], default: 0
% Do you want to plot only the identity per data set (0) or also the
% identity distribution per replicate (1) ?
ident_plot = 0;

% % Length histogram + exponential fit
% Do you want to exclude segments from the fit?
startbin_fit = 2; % default: 2 , i.e. exclude first bin 
max_lenfit = 16600; % [bp], default: 22000

% Would you like to exclude segments with an identity smaller min_id or
% larger max_id ?
min_id = 0;  % default: 0
max_id = 1; % default: 1
% Would you like to change the BinWidth of the length distribution hist?
binwidth_len = 1000; % default: 1000

% Percentage of genome replaced (= Transfer) fit 
% Do you want have the origin as fixed point (false) or leave the y-axis as 
% a free parameter (true)
intercept = true;

%% Loading and calculating ...

% % Read in Data, check if all you need is available and write into CNPplot

% Read in
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
    [log, pos_tmp] = ismember(dataset{i}, [Cluster(:).Sample]);
    ERR1 = find(log == 0);

    % Make a list of the sample positions in the Cluster struct (in the right
    % order)
    pos = [pos pos_tmp];

    if ~isempty(ERR1)
        error(['At least one sample of your data set ', num2str(i), ' is not in the loaded CNPSummaries.']);
    end
end
clear ERR1 pos_tmp log

% Collect only the clusters for this plot in CNPplot (and put them in 
% the right order)
CNPplot = Cluster(pos);
clear Cluster

% For better use, split the CNPsummaries - adist, ident
for i = 1 : length(CNPplot)
    if ~isempty(CNPplot(i).Adist)
        adist{i} = CNPplot(i).Adist(:,1);
    end
end
ident = {CNPplot(:).Ident};

% Preallocating ident_ds and adist_ds for speed
ident_ds = cell(1, numel(dataset)); 
adist_ds = cell(1, numel(dataset)); 

% Collect all identities and adists that belong to one data set
k = 1; % Counting position ins sampleP
for i = 1 : length(dataset)
    ident_ds{i} = vertcat(ident{k : k + length(dataset{i}) - 1});
    adist_ds{i} = vertcat(adist{k : k + length(dataset{i}) - 1});
    k = k + length(dataset{i});
end

%% Plotting ...

%
%   Figure 1 : Length distribution + Transfer per replicate
%
p1 = figure(1); hold on; 
set(gcf, 'renderer', 'painters', 'position', [990 490 length(sampleP)*150 500]);

for i = 1 : length(sampleP)
        subplot(1, length(sampleP), i); hold on; box on;
        set(gca, 'FontSize', 16);
        h_tmp(i) = histogram(adist{i}, 'BinWidth', binwidth_len, 'Orientation', 'horizontal', 'FaceColor', [0 0 0]);
        b_tmp_flip = barh(h_tmp(i).BinEdges(1:end-1) + h_tmp(i).BinWidth / 2, - h_tmp(i).Values , 'BarWidth', 1, 'FaceColor', [0.1 0.1 0.1]);
        set(gca, 'YTick', [], 'XTick', []);
        
        if i == 1
           ylabel('Length of integrated segments [bp]');
           set(gca, 'YTick', 0: 5E3: 60E3);
        end
        
        xlim([-max(h_tmp(i).Values)-1 max(h_tmp(i).Values)+1])
        ylim([0 round(max(vertcat(adist{:})), -3) + 2000]);
        xlabel(char(CNPplot(i).Sample))
        
        yyaxis right
        plot(0, CNPplot(i).Transfer, 'r*','MarkerSize', 10, 'LineWidth', 1.5)
        set(gca, 'YTick', [], 'XTick', [], 'YColor', 'r');
        if i == length(sampleP)
            if max([CNPplot(:).Transfer]) > 25
                steps = 5;
            elseif max([CNPplot(:).Transfer])> 10
                steps = 2;
            elseif max([CNPplot(:).Transfer]) > 5
                steps = 1;
            else 
                steps = 0.5;
            end
            set(gca, 'YTick', 0 : steps : 40 );
            ylabel('Transfer [%]');
        end
        ylim([0 round(max([CNPplot(:).Transfer]))+1])
end


%
%   Figure 2 : Identity plot - all replicates merged +/- single
%
% 

p2 = figure(2); hold on; box on;
title("All detected clusters > " + min_len + " bp")
if ident_plot == 1
    subplot(1,2,1); 
    hold on; box on;
    set(p2, 'renderer', 'painters', 'position', [10 10 1100 500],'paperpositionmode','auto'); 
elseif ident_plot == 0
    set(p2, 'renderer', 'painters', 'position', [10 10 700 500],'paperpositionmode','auto'); 
else
    error('To show the identity plot correctly, choose ident_plot = 1 or ident_plot = 0.');
end
set(gca, 'FontSize', 16);

% Plot 1: Identity histograms are shown per dataset
for i = 1 : length(dataset)
    % Create temporary variables ident_hist/adist_hists
    ident_hist = ident_ds{i};
    adist_hist = adist_ds{i};

    % Plot identity histogram of segments with minimum length min_len
    h(i) = histogram(ident_hist(adist_hist >= min_len)*100, 'BinWidth', 1, 'FaceColor', clr_ds(i,:), 'FaceAlpha', 0.2, 'EdgeColor', clr_ds(i,:), 'EdgeAlpha', 0.8, 'LineWidth', 1.5, 'Normalization', 'probability'); 

    % Calculate mean identity and the standard deviation
    mean_id(i) = mean( ident_hist(adist_hist >= min_len) * 100 );
    std_id(i) = std( ident_hist(adist_hist >= min_len) * 100 );
end

xlim([70 100]);
xlabel('Identity of integrated segments [%]');
ylim([0 max([h(:).Values])+0.01]);
ylabel('Probability');

leg = [];
for i = 1 : length(mean_id)
   leg{i} = [ds_name{i}, ': (', num2str(mean_id(i), '%0.1f'), ' \pm ', num2str(std_id(i), '%0.1f'), ') %'];
end

legend(char(leg), 'Location', 'nw');
clear leg

% Create a second plot with the identities for all samples seperately

if ident_plot == 1
    subplot(1,2,2); 
    hold on; box on; 
    set(gca, 'FontSize', 16);
    
    for i = 1 : length(CNPplot)
        tmp = ident{i};
        histogram(tmp(adist{i} >= min_len)*100, 'FaceColor', clr_samp(i,:), 'Normalization', 'probability', 'BinWidth', 1, 'FaceAlpha', 0.2, 'EdgeColor', clr_samp(i,:),'EdgeAlpha', 0.8, 'LineWidth', 1.5);
    end
        
    xlabel('Ident. of integr. segments per replicate [%]');
    ylabel('Probability');
    xlim([70 100]);
    legend(char(sampleP), 'Location', 'nw');
end


%
%   Figure 3: Length distributions + exponential fit
%

p3 = figure(3); hold on; box on;
set(p3, 'renderer', 'painters', 'position', [10 590 800 400],'paperpositionmode','auto'); 
set(gca, 'FontSize', 16);

for i = 1 : length(dataset)
    % Create a temporary variable for the current data set
    tmp = adist_ds{i};

    % Exclude the adists that have an identity higher than max_id or
    % smaller than min_id
    adist_hist = tmp(ident_ds{i} <= max_id & ident_ds{i} >= min_id);

    % Plot the adist histogram for the current dataset
    h3(i) = histogram(adist_hist, 'BinWidth', binwidth_len, 'FaceColor', clr_ds(i,:),'FaceAlpha', 0.2, 'EdgeAlpha', 0.8, 'EdgeColor', clr_ds(i,:), 'Normalization', 'probability', 'LineWidth', 1.5);
    
    leg3{i} = ds_name{i};
end
clear tmp

xlim([0 max(vertcat(adist{:})) + binwidth_len]);
xlabel('Length of integrated segments [bp]');
ylabel('Probability');
ylim([0 round(max([h3(:).Values]),2) + 0.02]);

%
%   Figure 32 : Lin. Fit
%
p32 = figure(32);
hold on; box on;
set(p32, 'renderer', 'painters', 'position', [600 10 800 500],'paperpositionmode','auto'); 
set(gca, 'FontSize', 16);
k = 1; %legend counter
for i = 1 : length(dataset)
    figure(32); hold on;
    y = h3(i).Values;
    y0 = y(y>0);
    logy0 = log(y0);
    x = h3(i).BinEdges(1:end-1)+h3(i).BinWidth/2;
    x0 = x(y>0);
    plot(x0,logy0, 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', clr_ds(i,:));
    cutidx = find(x0 > max_lenfit, 1, 'first');
    if isempty(cutidx); [~,cutidx] = max(x0); end
    ft32 = fitlm(x0(startbin_fit:cutidx),logy0(startbin_fit:cutidx), 'intercept', intercept)
    slope32(i) = ft32.Coefficients.Estimate(2);
    slope32err(i) = ft32.Coefficients.SE(2);
    interc32(i) = ft32.Coefficients.Estimate(1);
    xplot = 0:1000:round(max(x0),-3)+1000;
    yplot = xplot*slope32(i) + interc32(i);
    plot(xplot, yplot, 'LineWidth', 1.5, 'Color', [clr_ds(i,:) 0.5]);
    leg32{k} = ds_name{i};
    leg32{k+1} = ['Characteristic length = (', num2str(-log(2)/slope32(i), '%.0f'), ' \pm ', num2str(log(2)/slope32(i)^2*slope32err(i), '%.0f'), ')  bp'];
    k = k + 2;
    %
    %   Figure 33: Regression plot of figure(32) - fit 
    %
    figure(33); hold on; box on; grid on; 
    set(gcf, 'renderer', 'painters');
    xres = x0(x0<max_lenfit);
    yres = logy0(x0<max_lenfit) - (xres*slope32(i) + interc32(i));
    
    plot(xres, yres, 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', clr_ds(i,:) );
end
figure(33); hold on; set(gcf, 'renderer', 'painters');
xlabel('Length of integrated segments [bp]');
ylabel('Residual');
plot([0 max_lenfit+1000], [0 0], 'r--', 'LineWidth', 1.5);
legend(char(leg32([1:2:length(dataset)*2-1]))); % nimmt jeden zweiten Eintrag aus leg32

figure(32); hold on;
xlim([0 max([h3(:).BinEdges])]);
ylim([floor(log(min(nonzeros([h3(:).Values])))) 0]);
legend(char(leg32));
xlabel('Length of integrated segments [bp]');
ylabel('log(probability)');



%
% Plot the fit 32 in figure(3) 
%
figure(3); hold on; 
c = length(dataset) + 1; %legend counter
for i = 1 : length(dataset)
    yplot3 =exp(interc32(i))*exp(xplot*slope32(i));
    plot(xplot, yplot3, 'LineWidth', 2, 'Color', clr_ds(i,:))
    leg3{c} = ['Characteristic length = (', num2str(-log(2)/slope32(i), '%.0f'), ' \pm ', num2str(log(2)/slope32(i)^2*slope32err(i), '%.0f'), ')  bp'];
    c = c + 1;
end
legend(char(leg3), 'Location', 'ne');
clear leg3

%
% Figure 34: CDF + fit
%ecdf(y,'Bounds','on')
% hold on
% plot(x,evcdf(x,0,3))
% grid on
% title('Empirical CDF')
% legend('Empirical CDF','Lower Confidence Bound','Upper Confidence Bound','Theoretical CDF','Location','best')
% hold off

for i = 1 : size(dataset,2)
    figure(34+i-1); hold on; box on;
    set(gcf, 'renderer', 'painters', 'position', [1100 400 750 400],'paperpositionmode','auto'); 
    % Plot the empirical cdf (exp data)
    tmp_fd = adist_ds{i};
    cdfplot(tmp_fd)
    %ecdf(tmp_fd, 'Bounds', 'on')
    
    fd1 = fitdist(tmp_fd, 'exponential');
    xval = 1:20000;
    %plot(xval, evcdf(xval,0,3))
    ycdf = cdf(fd1, xval);
    plot(xval, ycdf, 'LineWidth', 1.5, 'Color', [0 1 0 0.4]);
    % Create a tmp for the fit data
    
    cdfplot(tmp_fd(tmp_fd>500)-500);
    %ecdf(tmp_fd(tmp_fd>500)-500, 'Bounds', 'on')
    fd2 = fitdist(tmp_fd(tmp_fd>500)-500, 'exponential');
    ycdf2 = cdf(fd2, xval);
    plot(xval, ycdf2, 'LineWidth', 1.5, 'Color', [1 0 0 0.5]);
    title(ds_name{i});
    xlim([0 2E4])
    set(gca, 'FontSize', 14);
    legend('CDF All Exp Data', ['L1/2 = ', num2str(log(2)*fd1.mu, '%0.0f'), ' bp'], 'CDF Exp Data > 500 bp', ['L1/2 = ', num2str(log(2)*fd2.mu, '%0.0f'), ' bp']);
end

%
%   Figure 4: Percentage of genome replaced  
%
% In this figure, we do not care about the predefined data sets. All
% samples are put into evol exp sets (each evolexp is a set - Ws, Wns, ...) 
% and the information about the cycle is extracted from the sample name

p4 = figure(4);
hold on; box on; grid on;
set(p4, 'renderer', 'painters', 'position', [1200 10 750 600],'paperpositionmode','auto'); 
set(gca, 'FontSize', 18);

% What variables do we already have for this plot ?
sampleP; % = name of the sample
transfer = [CNPplot(:).Transfer]; % = percentage of genome replaced

% Find out to which evol exp the samples belong, see variable evolexp
% Preallocating variables evolexp and set for speed
evolexp = cell(1, numel(sampleP));
setx = zeros(1, numel(sampleP));

for i = 1 : length(sampleP)
    if length(sampleP{i}) < 9
        evolexp{i} = sampleP{i}(1 : length(sampleP{i}) - 4);
    else
        evolexp{i} = sampleP{i}(1 : length(sampleP{i}) - 6);
    end
end

% How many different evol exps are in the datasets and what are their 
% prefixes?
evolexp_sets = unique(string(evolexp));
clr_ev = clrmp; % summer(length(evolexp_sets));

% Compare the samples to the found prefixes evolexp_sets
for j = 1 : length(sampleP)
    for i = 1 : length(evolexp_sets)
        if strcmp(evolexp_sets(i), evolexp{j})
            setx(j) = i;
        end
    end
end

% And after how many cycles was this sample sequenced?
% Preallocate variable cy for speed
cy = zeros(1, numel(sampleP));

for i = 1 : length(sampleP)
    % Extract the cycle number from the sample name and wrtie it in var. cy
    if length(sampleP{i})<9
            cy(i) = str2double(sampleP{i}(end-1:end));  
    else
            cy(i) = str2double(sampleP{i}(end-3:end-2));  
    end
end

% And now, plot !
% The following part is quite tricky, because it fits for every evol exp
% separately linear regressions to all possible combinations of your
% data, no matter how many cycles you included. I hope the comments help to
% make it easier for you to follow. 
% transfer == % of genome replaced

% Look at every evol exp (set) 'i' separately
for i = 1 : length(evolexp_sets)
    
    % Define temporary variables tmp_x and tmp_y for your evol exp i,
    % containing the percentage of genome replaced (tmp_y) and the
    % belonging cycl(tmp_x):
    tmp_x = cy(setx==i);
    tmp_y = transfer(setx==i);
    
    % Plot the data point of the current evol exp
    p44(i*2 - 1) = plot(tmp_x*2, tmp_y, 'x', 'MarkerEdgeColor', clr_ev(i,:), 'LineWidth', 1.5, 'MarkerSize', 8);
    
    % Now, we start with the preparations for the fits. We need to count
    % the number of possible combinations of transfer data, i.e. how many fits
    % do we need to make, and we create a matrix, with all the combinations
    % 'variants'.
    
    % What different cycles are in the data sets ? - cy_idx
    cy_idx = unique(tmp_x);
    
    % In the following loop, the tranfer values per cycle are counted
    % (length(tmp_trans)) and the number of possible transfer variants = 
    % the number of fits we are gonna make are counted.
    variants_no = 1;
    for j = 1 : length(cy_idx)
        % Write all transfer values (% of genome replaced) in a cell 
        % according to their cycle
        tmp_trans{j} = tmp_y(tmp_x == cy_idx(j));    
        % Calculate the number of possible variants (for the fit):
        variants_no = variants_no * length(tmp_trans{j});
    end
    
    % We create a temporary cy variable for the while-loop, every cycle
    % that is already in the variants matrix is deleted from the tmp_cy_idx
    tmp_cy_idx = cy_idx;

    % This variable count the number of tranfer values that are written in
    % the variants matrix.
    sum_len = 1;
    
    % Preallocate the variants matrix :
    variants = zeros(variants_no, length(cy_idx));
    
    % Now, recursively go through all different cycles, and create the
    % variants matrix
    while ~isempty(tmp_cy_idx)
        
        % current cycle
        cy_tmp = tmp_cy_idx(end);
        
        % All values belonging to this cycle
        val = tmp_y(tmp_x == cy_tmp);
        
        % Number of entries for this cycle
        val_len = length(val);
        
        % Where are we in the variants matrix right now ? - lincount
        lincount = 1;
        
        % As long as lincount is smaller than length of the variants matrix
        while lincount < length(variants)
            
            % We fill the last empty row with the transfer values of this
            % cycle. In the first iteration, sum_len equals 1, so very
            % value is written in one row and then it start again from the
            % beginning. In the second iteration, sum_len is 1 * # transfer
            % values of the cycle before, so the first transfer value is
            % written in sum_len rows and then the second as well until the
            % end of the variants matrix ..and so on. Like this, all
            % possible combinations, no matter how many cycles you give me,
            % can be found in the variants matrix.
            for k = 1 : length(val)
                tmp_var = ones(sum_len, 1) * val(k);
                variants(lincount : lincount + sum_len - 1, length(tmp_cy_idx)) = tmp_var;
                lincount = lincount + sum_len;
            end
        end
        
        % As described in the paragraph above, the sum_len is increased.
        sum_len = sum_len * length(val);
        
        % And the last entry of the temporary cycles array is deleted.
        tmp_cy_idx(end) = [];
        
        % The loop stops, when tmp_cy_idx is empty.
    end    
   
    % Now we finally can fit! 
    % We go through all variants and do a linear regression
    % First add zero as fit point
    variants = [variants zeros(length(variants),1)];
    for n = 1 : length(variants)
        ft = fitlm([cy_idx*2 0], variants(n,:), 'intercept', intercept);
        
        % When you decided to use the origin as fixed point, intercept ==
        % 0, otherwise , the intercept is written into interc - array:
        if intercept == 1
            slope(n) = ft.Coefficients.Estimate(2);
            interc(n) = ft.Coefficients.Estimate(1);
        elseif intercept == 0
            slope(n) = ft.Coefficients.Estimate;
            interc = 0;
        end
    end
    
    % And plot
    figure(4); hold on; box on;
    % Calculate the y value from the mean of slope (and intercept) and the std    
    xmean = [0 max(cy)*2 + min(cy)*2];

            ymean = mean(slope) * xmean + mean(interc);
            ystdpos = (mean(slope) + std(slope)) * xmean + mean(interc) + std(interc);
            ystdneg =(mean(slope) - std(slope)) * xmean + mean(interc) - std(interc);
        
        p44(i*2) = plot(xmean, ymean, '-', 'color', [clr_ev(i,:) 0.8], 'LineWidth', 1.5);
        p42(i) = plot(xmean, ystdpos, '--', 'color',[clr_ev(i,:) 0.4], 'LineWidth', 1.5);
        p43(i) = plot(xmean, ystdneg, '--', 'color',[clr_ev(i,:) 0.4], 'LineWidth', 1.5);

        
        leg{i*2 - 1} =  char(evolexp_sets(i));
        leg{i*2} = ['(',num2str(mean(slope),'%.2f'),'\pm', num2str(std(slope),'%.2f'),') % h^{-1}'];

end

l4 = legend(p44(:), char(leg), 'Location', 'nw', 'LineWidth', 1.2);
xlabel('DNA uptake [h]');
ylabel('Genome replaced [%]');
xlim([0 max(cy)*2 + min(cy)*2]);
ylim([0 max(transfer) + min(transfer)]);

% Fig 1 to the front
figure(1) 
close 1
figure(2)
%% Saving ...

yes = 1; no = 0; 
% saveme = input('Would you like to save the plots in the savepath ? (yes/no) ');
saveme = no;
if isempty(saveme) || saveme == 0
    disp('Ok, the plots are not saved.');
elseif saveme == 1
    %print(p1, '-painters', '-dpng', [savepath,datestr(now, 'yyyymmdd'),'_', 'TransferPlot_', [sampleP{:}]]);
    print(p2, '-painters', '-dpng', [savepath,datestr(now, 'yyyymmdd'),'_', 'IdentHist_', char(evolexp_sets)]); 
    print(p3, '-painters', '-dpng', [savepath,datestr(now, 'yyyymmdd'),'_', 'AdistHist_',  char(evolexp_sets)]);
    print(p32,  '-painters', '-dpng', [savepath,datestr(now, 'yyyymmdd'),'_', 'AdistHistLog_',  char(evolexp_sets)]);
    print(p4, '-painters', '-dpng', [savepath,datestr(now, 'yyyymmdd'),'_', 'TransferFit_',  char(evolexp_sets)]);
    % print(figure(34), '-painters', '-dpng', [savepath,datestr(now, 'yyyymmdd'),'_', 'CDF_LenDistr_',  char(evolexp_sets), '_cy10']);
    % print(figure(35), '-painters', '-dpng', [savepath,datestr(now, 'yyyymmdd'),'_', 'CDF_LenDistr_',  char(evolexp_sets), '_cy20']);
    disp('All plots have been saved.');
else
    disp('Ok, the plots are not saved.');
end