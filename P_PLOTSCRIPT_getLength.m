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
    time_mask = ~contains([Data.CNPSummary.Sample],"10");
    Cluster = [Cluster; Data.CNPSummary(time_mask)];
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