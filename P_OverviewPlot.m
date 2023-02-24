%
%
%%%%%%%%%%%%%%%%   OverviewPlot created by MonaIsa  %%%%%%%%%%%%%%%%%%%%%%%%
%
%  What this script does: It plots, along the whole genome, the different 
%      genomic changes that have occurred such as recombinations, denovos,  
%      deletions and insertions. The user can 
%
%  Input: Cluster information; operon list; bed file (; masterlist)
%    -> the user can load a number of SNP2CNP outputs, they are collected
%    automatically and the program, when plotting, merely searches for the 
%    sample names
%
%  Output: OverviewPlot

clear all; close all

%% give the paths/name to the acc genome, the data and give the samples you want to plot
%Here you can add as many Cluster information as you want
basePath = "/home/isabel/Documents/Doktorarbeit/P_Populations/summaries/";

CNPname{1} = basePath + "20220506_Pop_LibBvalwS2_CNPSummary.mat";
CNPname{2} = basePath + "20220506_Pop_LibBvalwS3_CNPSummary.mat";
CNPname{3} = basePath + "20220509_LibBvalwS2d12_CNPSummary_withCluFreq.mat";
CNPname{4} = basePath + "20220509_LibBvalwS3d5_CNPSummary_withCluFreq.mat";

%Here add the deldup Outputs.
DDname{1} = "/home/isabel/sciebo/ResultsShared/DFE_HighThroughput/2022_withSelection/DNAseq/20220425_LibBvalwS2d12_deldup.mat";
%DDname{2} = basePath + "...mat";

%Acc genome
acc = '/home/isabel/sciebo/ResultsShared/kleinesPaper/allLists/accBval2Bs166NCe.txt';

% what are the samples you want to plot in this order from top to bottom?
% (the names must be the same as are given in the SNP2CNP output structure
sampleP = ["LibBval_Pop_S3d3" "LibBval_Pop_S3d5" "LibBvalwS3d5_56"]

%User input, in which row from the top would you like to plot which
%feature?
clu    = [1 1 1];
del    = [0 0 0];
dup    = [0 0 0];
denovo = [0 0 0];

% How should the plot look like?
recipsize = 4215607;
rows = 3;
region = [3555000 3700000];
if isempty(region)
    region = [1 recipsize];
end

% If wanted: sample name in the plot, if different from sampleP
nameinP = ["PopBVALevoCM d3" "PopBVALevoCM d5" "BVALevoCM"];

SavePlot = "OFF";
PlotName = "BVALevoDMPop";
savepath = "/home/isabel/Documents/Doktorarbeit/P_Populations/plots/";

%% load the data
% load the data that will be plotted
Cluster = [];
for i=1:size(CNPname,2)
    Data = load(CNPname{i});
%     if i==1 || i== 2
%         Data.CNPSummary = rmfield(Data.CNPSummary,{'ClusterFrequency','Info'});
%     end
    Cluster = [Cluster;Data.CNPSummary];
end
clear Data
deletions = []; duplications = [];
for i=1:size(DDname,2)
    Data = load(DDname{i});    
    deletions = [deletions;Data.deldup.del'];
    duplications = [duplications;Data.deldup.dup'];
end
clear Data
%Accessory genome
fid = fopen(acc);
refchr = textscan(fid,'%f %f %f','delimiter','\t');
fclose(fid);
Acc.S=refchr{1}; Acc.E=refchr{2};
clear refchr

%% Initialize the plot and add Acc. genome
f12=figure('Units','centimeters','position',[20 1 16.97 24],'paperpositionmode','auto');   
    for j = 1:size(Acc.S,1)
        for b = 1:rows
            subplot(rows,1,b)
            plot([Acc.S(j) Acc.E(j)],[0 0],'Color',[0.7 0.7 0.7],'linewidth',6)
            hold on
        end
    end
    
samnum = numel(sampleP);
start = region(1);
edge = region(2);
ell = (edge-start+1)/(rows - (1/8));
    for b = 1:rows
        subplot(rows,1,b)
        xlim([start+ell*(b-1),start+ell*b]); 
        
        ylim([-0.5 samnum+0.5]); ax = gca;
        ax.YTick=0:(samnum+1);
        set(gca,'TickLabelInterpreter','none')
        if isempty(nameinP)
            names = flip(sampleP)
            ax.YTickLabel={'B168 Acc' names{:}};
        else
            names = flip(nameinP)
            ax.YTickLabel={'Acc genome' names{:}};
        end
        ax.XGrid = 'off';
        ax.YGrid = 'on';
    end
ax = gca; ax.XAxis.Exponent=6; clear b j ax i 
subplot(rows,1,rows)
xlabel('Chromosome position with respect to \it{B. subtilis}');
xlim([start+ell*(rows-1),edge]);

%% Recombinations are plotted, color matches length

% Here, the maximum length in all samples is searched to set a color limit
if sum(clu~=0)>0
    Cdistcollect = []; VarFreq = [];
    for i=1:numel(clu)
        if clu(i)~=0
            pos=find(ismember([Cluster(:).Sample],sampleP(i)));
            Cdistcollect = [Cdistcollect;Cluster(pos).Cdist];
%             VarFreq = [VarFreq;Cluster(pos).];
        end
    end
%set colors
    clrlimit = max(Cdistcollect);
    clrmap=(jet(clrlimit));
    clrmap(end+1,:)=[0 0 0];
    
 maskall = zeros(recipsize,1); 
    for i = 1:numel(sampleP)
        if clu(i)~=0
            L = samnum - i +1;
            clear idx CS CE CL
            idx = find(ismember([Cluster(:).Sample],sampleP(i)));
            if isempty(Cluster(idx).C)
                continue
            end
            CS = Cluster(idx).C(1,:);
            CE = Cluster(idx).C(2,:);
            CL = Cluster(idx).Cdist;
            mask = zeros(recipsize,1);
             
            figure(f12); hold on;
            for j = 1:size(CS,2)
                if CS(j)>CE(j)
                    for b = 1:rows
                        mask(round(CS(j)):recipsize) = ones(recipsize-round(CS(j))+1,1);
                        mask(1:round(CE(j))) = ones(round(CE(j)),1);
                        subplot(rows,1,b)
                        hold on;
                        plot([CS(j) recipsize],[L L],'color',clrmap(CL(j),:),'linewidth',6)
                        plot([1 CE(j)],[L L],'color',clrmap(CL(j),:),'linewidth',6)
                    end
                else
                    for b = 1:rows
                        mask(round(CS(j)):round(CE(j))) = ones(round(CE(j))-round(CS(j))+1,1);
                        subplot(rows,1,b)
                        hold on;
                        plot([CS(j) CE(j)],[L L],'color',clrmap(CL(j),:),'linewidth',6)
                    end
                end
            end
            maskall = maskall + mask;
        end       
    end
stepsize = clrlimit/2;    
    subplot(rows,1,rows)
    ticks = 0:stepsize:clrlimit;
    tickLabel = cellstr(string(ticks/1000));
    colorbar('eastoutside','TickLabels',{tickLabel{1:end-1} strcat(tickLabel{end},'x10^3') 'interpreter' 'latex'}, ...
        'Ticks', 0:stepsize:clrlimit);
    cmap = colormap(jet);
    caxis([0, clrlimit]);
end
%% Denovos
if sum(denovo~=0)>0 
    for i = 1:numel(sampleP)
        L = samnum - i +1;
        idx = find(ismember([Cluster(:).Sample],sampleP(i)));
        if denovo(i)==1
            deno = [];
            deno = Cluster(idx).denovo;
            for j=1:size(Cluster(idx).denovo,2)
                for b = 1:rows
                    subplot(rows,1,b)
                    if deno(1,j)>=start+ell*(b-1) && deno(1,j)<=start+ell*b && deno(1,j)<= edge
                    hold on;                    
                    text(deno(1,j),L,'|')
                    end
                end
            end
        end
    end
end

%% Deletions and Duplications

if sum(del~=0)>0
    for i = 1:numel(sampleP)
        L = samnum - i +1;
        
        idx = find(contains({deletions(:).sample},sampleP{i}));
        if del(i)~=0
            de_start = []; de_edge = [];
            de_start = deletions(idx).start;
            de_edge = deletions(idx).edge;
            for j=1:size(de_start,2)
                for b = 1:rows
                    subplot(rows,1,b)
                    plot([de_start;de_edge],[L*ones(numel(de_start),1)';L*ones(numel(de_start),1)'],'k','LineWidth',6)
%                     if de_start(1,j)>=start+ell*(b-1) && de_start(1,j)<=start+ell*b && de_start(1,j)<= edge                       
%                         text(de_start(1,j),L,'[','Color','k','FontSize',14)
%                     end
%                     if de_edge(1,j)>=start+ell*(b-1) && de_edge(1,j)<=start+ell*b && de_edge(1,j)<= edge                       
%                         text(de_edge(1,j),L,']','Color','k','FontSize',14)
%                         hold on
%                     end
                end
            end
        end
    end
end

if sum(dup~=0)>0
    for i = 1:numel(sampleP)
        L = samnum - i +1;
        idx = find(contains({duplications(:).sample},sampleP{i}));
        if dup(i)~=0
            duo_start = []; dup_edge = [];
            if isempty(duplications(idx).start)
                continue
            end
            dup_start = duplications(idx).start;
            dup_edge = duplications(idx).edge;
            for j=1:size(dup_start,2)
                for b = 1:rows
                    subplot(rows,1,b)
                    plot([dup_start;dup_edge],[L*ones(numel(dup_start),1)';L*ones(numel(dup_start),1)'],'m','LineWidth',6)
                    if dup_start(1,j)>=start+ell*(b-1) && dup_start(1,j)<=start+ell*b && dup_start(1,j)<= edge                       
                        text(dup_start(1,j),L,'{','Color','m','FontSize',14)
                    end
                    if dup_edge(1,j)>=start+ell*(b-1) && dup_edge(1,j)<=start+ell*b && dup_edge(1,j)<= edge                       
                        text(dup_edge(1,j),L,'}','Color','m','FontSize',14)
                        hold on
                    end
                end
            end
        end
    end
end


%% HotCold per basepair - make overview plot for clusters
figure(4)
subplot(6,1,1)
rows = 6; 
pstart = region(1); pedge = region(2); plength = pedge-pstart+1; step = plength/rows;
title('Hot and Cold positions per bp')
hold on
for i=1:rows
    subplot(rows,1,i)
    %reg = [floor(1+(i-1)*(pos/rows)):ceil((pos/rows)*i)];
    reg = [floor(pstart+(i-1)*step):ceil(pstart+step*i-1)];
    bar(reg,maskall(reg),'BarWidth',1)
    %xlim([1+(i-1)*(pos/rows) (pos/rows)*i ])
    xlim([pstart+(i-1)*step pstart+step*i])
    ylim([0 max(maskall)+1])
end
subplot(rows,1,ceil(rows/2))
ylabel('# Hits')
subplot(rows,1,rows)
xlabel('Position on Genome');

%% print plot


if 	SavePlot == "ON"
    print(gcf, '-painters', '-dpng', [savepath + datestr(now, "yyyymmdd") + "_CM_OverviewPlot_" + PlotName]);
else
    disp('Ok, the plot is not saved.')
end
