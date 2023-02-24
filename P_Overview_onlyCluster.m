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

basePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/DNASeq/2a_Cluster/";
CNPname{1} = basePath + "20220630_Vns_CNPSummary.mat";


%Acc genome
acc = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/acc/accBval_2NCe.mat";

% what are the samples you want to plot in this order from top to bottom?
% (the names must be the same as are given in the SNP2CNP output structure
%sampleP = ["LibBvalwS2d12_6","LibBvalwS2d12_14","LibBvalwS2d12_16","LibBvalwS2d12_30","LibBvalwS2d12_38","LibBvalwS2d12_39"];
sampleP = ["Vns0410" "Vns0810" "Vns0420" "Vns0820.5"];

%User input, in which row from the top would you like to plot which
%feature?
clu    = [1 1 1 1];
denovo = [0 0];

% How should the plot look like?
recipsize = 4215607;
rows = 3;
region = [];
if isempty(region)
    region = [1 recipsize];
end

% If wanted: sample name in the plot, if different from sampleP
nameinP = [];

SavePlot = "OFF";
savepath = "/media/isabel/IsabelFestplatte/Doktor/MATLAB_plots/";

%% load the data
% load the data that will be plotted
Cluster = [];
for i=1:size(CNPname,2)
    Data = load(CNPname{i});
    Cluster = [Cluster;Data.CNPSummary];
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
        xlim([start+ell*(b-1),start+ell*b]); grid on;
        ylim([-0.5 samnum+0.5]); ax = gca;
        ax.YTick=0:(samnum+1);
        set(gca,'TickLabelInterpreter','none')
        if isempty(nameinP)
            names = flip(sampleP);
            ax.YTickLabel=["Bs166 Acc" names];
        else
            names = flip(nameinP);
            ax.YTickLabel=["Acc genome" names];
        end
    end
ax = gca; ax.XAxis.Exponent=6; clear b j ax i 
subplot(rows,1,rows)
xlabel('Genome position with respect to \it{B. subtilis}');
xlim([start+ell*(rows-1),edge]);

%% Recombinations are plotted, color matches length

% Here, the maximum length in all samples is searched to set a color limit
if sum(clu~=0)>0
    Cdistcollect = [];
    for i=1:numel(clu)
        if clu(i)~=0
            pos=find(ismember([Cluster(:).Sample],sampleP(i)));
            Cdistcollect = [Cdistcollect;Cluster(pos).Cdist];
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
                        plot([CS(j) recipsize],[L L],'color',clrmap(CL(j),:),'linewidth',4)
                        plot([1 CE(j)],[L L],'color',clrmap(CL(j),:),'linewidth',4)
                    end
                else
                    for b = 1:rows
                        mask(round(CS(j)):round(CE(j))) = ones(round(CE(j))-round(CS(j))+1,1);
                        subplot(rows,1,b)
                        hold on;
                        plot([CS(j) CE(j)],[L L],'color',clrmap(CL(j),:),'linewidth',4)
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
        if denovo(i)~=0
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
    reg = round([floor(pstart+(i-1)*step):ceil(pstart+step*i-1)]);
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

if SavePlot == "ON"
    print(gcf, '-painters', '-dpdf', [savepath, 'OverviewPlot_', [names{:}]]);
else
    disp('Ok, the plot is not saved.')
end
