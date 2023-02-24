%
%
%%%%%%%%%%%%%%%%   Cov2DelDup created by Isabel  %%%%%%%%%%%%%%%%%%%%%%%%
%
%
% What this script does: It finds all the deletions and duplications for
% given coverage data
%
%    input:
% -- is the coverage of the samples of interest; the artefacts to be excluded
% -- and the parameters of choice
%
%    output:
% -- is a structure deletion and duplication - together saved as deldup
% -- this script offers to look at overview plots that can help you to
%    decide if the positions you found are reasonable

clear all; close all
%% give the paths/names and set the parameters
% ---------------------------load the coverage

% covpath = "/home/isabel/sciebo/ResultsShared/DFE_HighThroughput/DNASeq_LibBmoj/";
% covsample = "LibSCBmoj" + string(1:10);
covpath   = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper_BigData/Hybrids2Bs166/coverage/";
%covpath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper_BigData/Bs1662Donor/";
%covsample = "Bs166";
covsample =  "Geons0" + ["110" "210" "310" "410" "510" "610" "710" "810" "120" "220" "320" "420" "520"];
%covsample = "BAns0" + ["110" "210" "310" "410" "610" "710" "810" "120" "220" "320" "420" "620" "720" "820"]; % 
%covsample = "Vns0" + ["110" "210" "310" "410" "510" "610" "710" "810"...
%    "120" "220.5" "320" "420" "520.5" "620.5" "720.5" "820.5"]; % [1:5]
%covsample = "Geons" + ["0620.5" "0720.5" "0820.5"];
covsuffix = "_2NCe_coverage.txt";

expName = "Geons";


clrmap = jet(numel(covsample));

% variables 
refchr=4215607; 
cutoff = 5;       % this is the factor c: cutoff= mu-c*sigma
sldw = 30;        % sliding window that smooths the cov prior to analysis
minL = 1;         % minimal length of detected segment
del_para = 0;     % this is the absolute value that del has to hit
min_consdel = 10; % del_para has to be hit this number of consecutive times
dup_para = 2; % this is the factor to the mean that duplication has to hit
min_consdup = 10; % this is how often dub_para has to be hit consecutively

% SET PREFERENCES!
exclude     = "ON";
saveOut     = "ON";
plotSamples = "ON";

%%
% load artefact regions and exclude them
arte.path = '/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/Bs166NCe_ArteCov.txt';
if ~isempty(arte.path)
fid = fopen(arte.path);
imp = textscan(fid,'%f %f %s %s','headerLines', 1);
fclose(fid);
arte.InS = imp{1};
arte.InE = imp{2};
arte.type = imp{3};
end
%load mm regions and exclude them
mm.path = '/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/Bs166NCe_mm.txt';
if ~isempty(mm.path)
fid = fopen(mm.path);
imp = textscan(fid,'%f %f %f');
fclose(fid);
mm.InS = imp{1};
mm.InE = imp{2};
end



%% preparations: 
% load the data
for i=1:numel(covsample)
    fid = fopen([covpath + covsample(i) + covsuffix]);
    imp = textscan(fid,'%s %f %f');
    fclose(fid);
    in(i).cov = imp{3};
    in(i).sample = covsample{i};
    in(i).mu = mean(in(i).cov);
    in(i).sigma = std( in(i).cov);
    in(i).pos = 1:numel(in(i).cov);
    
    %plot the histograms as an overview
    figure(1)
    h(i)=histogram(in(i).cov,'FaceColor',clrmap(i,:),'EdgeColor',clrmap(i,:),'FaceAlpha',0.2,'BinWidth',5)
    hold on
    xlim([0 1000]) 
end
legend(h,covsample)
title('Overview of the coverage distributions')

% make an artefacts mask for del and dup
artedel_mask = ones(refchr,1); artedup_mask = ones(refchr,1);
if ~isempty(arte.path)
delidx = find(ismember({arte.type{:}},'del'));
dupidx = find(ismember({arte.type{:}},'dup'));
end
if ~isempty(mm.path)
edel = [arte.InS(delidx) arte.InE(delidx); mm.InS mm.InE];
edup = [arte.InS(dupidx) arte.InE(dupidx); mm.InS mm.InE];
end
if ~isempty(arte.path) || ~isempty(mm.path)
for i=1:size(edel,1)    
    artedel_mask(edel(i,1):edel(i,2)) = zeros(edel(i,2)-edel(i,1)+1,1);
end
for i=1:size(edup,1)  
    artedup_mask(edup(i,1):edup(i,2)) = zeros(edup(i,2)-edup(i,1)+1,1);
end
end
clear delidx dupidx

%% find the deletions
m = 0; 
for i=1:numel(covsample)
    m = m+1;
    %first set the paramters that are important and save them
    deletion(m).mu = mean(in(i).cov);
    deletion(m).sigma = std(in(i).cov);
    deletion(m).cutoff_factor = cutoff;
    if deletion(m).mu - cutoff*deletion(m).sigma >= 50
        deletion(i).thresh = deletion(m).mu - cutoff*deletion(m).sigma;
    elseif deletion(m).mu - cutoff*deletion(m).sigma < 50
        deletion(i).thresh = deletion(m).mu / 2;
        deletion(i).exception = true;
    end
    deletion(m).minL = minL;
    deletion(m).sldw = sldw;
    deletion(m).delpara = del_para; 
    deletion(m).minconsdel = min_consdel;
    % create the deletion mask and find deletions
    delmask0 = movmean(in(i).cov,sldw) <= deletion(i).thresh; % 1:deletions
    if exclude == "ON"
       delmask = double(delmask0) .* artedel_mask;  
    else 
        delmask = delmask0;
    end
    % Find the start of every deletion by shifting the list and adding it ontop
    sta = delmask(1:end) + [delmask(end);delmask(1:end-1)];
    sta2 = sta./[delmask(end);delmask(1:end-1)];
    start = in(i).pos(isinf(sta2));
    % Finding the ends of the deletions
    ed = delmask(1:end) + [delmask(2:end);delmask(1)];
    ed2 = ed./[delmask(2:end);delmask(1)];
    edge = in(i).pos(isinf(ed2));    
    % is there an ORI crossing?
    ORI = false;
    if ~isempty(edge) && edge(1) < start(1)
        ORI = true;
        deletion(m).ORI = 1;
        edge = [edge(2:end) edge(1)];
    end    
    % cut the deletions that are smaller than a number minL
    maskL = (edge - start +1) > minL;
    if ORI==1 && (edge(end) + (refchr - start(end) +1)) > minL
    maskL(end) = true;
    end
    start_minL = start(maskL);
    edge_minL = edge(maskL);   
    % now make sure that deletions contain at least min0 cons. 0
    % first make a mask with all 0 positions (they are 1)
    delmask2 = (movmean(in(i).cov,sldw) <= deletion(m).delpara);  
    stapel = delmask2;
    for p=1:min_consdel-1
        stapel =  stapel + [delmask2(end-p+1:end);delmask2(1:end-p)];
    end
    % go through stapel and check if min_consdel is fulfilled
    maskmin0 = false(numel(start_minL),1);
    for j=1:numel(start_minL)
        if ~isempty(find(stapel(start_minL(j):edge_minL(j)),10))
            maskmin0(j) = true;
        elseif ORI==1 && j==numel(start_minL) && (~isempty(find(stapel(start_minL(j):refchr),10)) || ...
                ~isempty(find(stapel(1:edge_minL(j)),10)))
            maskmin0(j) = true;
        end
    end
    deletion(m).sample = covsample(i);
    deletion(m).start = start_minL(maskmin0);
    deletion(m).edge = edge_minL(maskmin0);    
    % save other information 
    deletion(m).numD = size(deletion(m).start,2);% - mockdel.numD;
    deletion(m).sumL = sum(deletion(m).edge - deletion(m).start +1);% - mockdel.refLsum;
end
%clean up if wanted
clear an delmask delmask2 edge edge_minL fid imp j m maskL maskmin0 ORI p stapel start start_minL sta sta2 start ed ed2 edge

%% find the duplicatons
m = 0;
for i=1:numel(covsample)
    m = m+1;
     %first set the paramters that are important and save them
    duplication(m).mu = mean(in(i).cov);
    duplication(m).sigma = std( in(i).cov);
    duplication(m).cutoff_factor = cutoff;   
     if duplication(m).mu + cutoff * duplication(m).sigma <= duplication(m).mu*2 - 50
        duplication(i).thresh = duplication(m).mu + cutoff * duplication(m).sigma;
    elseif duplication(m).mu + cutoff * duplication(m).sigma > duplication(m).mu*2 - 50
        duplication(i).thresh = duplication(m).mu * 1.5;
        duplication(i).exception = true;
     end   
    duplication(m).minL = minL;
    duplication(m).sldw = sldw;
    duplication(m).duppara = dup_para * duplication(m).mu;
    duplication(m).minconsdup = min_consdup;
    % now create the duplication mask and find the dups
    dupmask0 = movmean(in(i).cov,sldw) >= duplication(i).thresh; % 1: duplications
    if  exclude == "ON"
       dupmask = double(dupmask0) .* artedup_mask;     
    else 
        dupmask = dupmask0;
    end
    % Find the start of every duplication by shifting the list and adding it ontop
    sta = dupmask(1:end) + [dupmask(end);dupmask(1:end-1)];
    sta2 = sta./[dupmask(end);dupmask(1:end-1)];
    start = in(i).pos(isinf(sta2));
    % Finding the ends of the duplications
    ed = dupmask(1:end) + [dupmask(2:end);dupmask(1)];
    ed2 = ed./[dupmask(2:end);dupmask(1)];
    edge = in(i).pos(isinf(ed2));
    % is there an ORI crossing?
    ORI = false;
    if ~isempty(edge) && edge(1) < start(1)
        ORI = true;
        duplication(m).ORI = 1;
        edge = [edge(2:end) edge(1)];
    end    
    % cut the duplications that are smaller than a number minL
    maskL = (edge - start +1) > minL;
    if ORI==1 && (edge(end) + (refchr - start(end) +1)) > minL
    maskL(end) = true;
    end
    start_minL = start(maskL);
    edge_minL = edge(maskL);   
    % now make sure that duplications contain at least min2mdub cons.
    % double values of the mean
    % first make a mask with all 0 positions (they are 1)
    dupmask2 = (movmean(in(i).cov,sldw) >= duplication(m).duppara);   
    stapel = dupmask2;
    for p=1:min_consdup-1
        stapel =  stapel + [dupmask2(end-p+1:end);dupmask2(1:end-p)];
    end
    % go through stapel and check if in the found deletions there is immer
    % at least once == min0
    maskmin0 = false(numel(start_minL),1);
    for j=1:numel(start_minL)
        if ~isempty(find(stapel(start_minL(j):edge_minL(j)),10))
            maskmin0(j) = true;
        elseif ORI==1 && j==numel(start_minL) && (~isempty(find(stapel(start_minL(j):refchr),10)) || ...
                ~isempty(find(stapel(1:edge_minL(j)),10)))
            maskmin0(j) = true;
        end
    end
    duplication(m).sample = covsample(i);
    duplication(m).start = start_minL(maskmin0);
    duplication(m).edge = edge_minL(maskmin0);    
    % save other information
    duplication(m).numD = size(duplication(m).start,2);% - mockdel.numD;
    duplication(m).sumL = sum(duplication(m).edge - duplication(m).start +1);% - mockdel.refLsum;
end
%clean up if wanted
clear an dupmask dupmask2 edge edge_minL fid imp j m maskL maskmin0 ORI p stapel start start_minL sta sta2 start ed ed2 edge

%% save data
if saveOut == "ON"

    deldup.del = deletion;
    deldup.dup = duplication;
    save(covpath + datestr(now,'yyyymmdd') + "_DelDup_" + expName + '.mat','deldup');
else
    disp('Ok, deletions and duplications are not saved.')
end

%% HotSpots for DelDup
delsum = zeros(refchr,1);
for i=1:numel([deletion(:).thresh])
    maskdels = zeros(refchr,1);
    for j=1:numel(deletion(i).start)
        maskdels(deletion(i).start(j):deletion(i).edge(j)) = 1;
    end
    delsum = delsum + maskdels;
end
dupsum = zeros(refchr,1);
for i=1:numel([duplication(:).thresh])
    maskdups = zeros(refchr,1);
    for j=1:numel(duplication(i).start)
        maskdups(duplication(i).start(j):duplication(i).edge(j)) = 1;
    end
    dupsum = dupsum + maskdups;
end
% all start and end points collected
delallstart = [deletion(:).start]; 
dupallstart = [duplication(:).start]; 

% plot the HotSpot plot
figure(33)
rows = 5;
pstart = 1;
pedge = refchr;
plength = pedge-pstart+1;
subplot(6,1,1)
for i=1:rows
    subplot(rows,1,i)
    region = ceil([pstart+(i-1)*(plength/rows):pstart+(plength/rows)*i-1]);
    bar(region,delsum(region),'b','HandleVisibility', 'off', 'BarWidth', 1)
    hold on
    bar(region,dupsum(region),'m','HandleVisibility', 'off', 'BarWidth', 1)
    if i==1
        plot(delallstart,max(delsum)*ones(numel(delallstart),1)+1,'Marker','d','MarkerFaceColor','b', 'LineStyle', 'none','Displayname','Deletion')
        plot(dupallstart,max(delsum)*ones(numel(dupallstart),1)+1,'Marker','*','MarkerFaceColor','m', 'MarkerEdgeColor','m','LineStyle','none','Displayname','Duplication')
    else
        plot(delallstart,max(delsum)*ones(numel(delallstart))+1,'Marker','d','MarkerFaceColor','b', 'LineStyle', 'none','HandleVisibility','off' )
        plot(dupallstart,max(delsum)*ones(numel(dupallstart))+1,'Marker','*','MarkerFaceColor','m', 'MarkerEdgeColor','m','LineStyle', 'none','HandleVisibility','off' )
    end
    xlim([pstart+(i-1)*(plength/rows) pstart+(plength/rows)*i])
    set(gca, 'XMinorTick', 'on', 'TickDir', 'out', 'TickLength', [0.005 0.025]);
    ylim([0 max(delsum)+2])
end
subplot(rows,1,1)
title('Hot and Cold plot')
legend('Location','NorthEastoutside')

%% Overview Plot
if plotSamples == "ON"
    rows = 5;
    pstart = 0;
    pedge = refchr;
    plength = pedge-pstart+1;
    for s=1:numel([deletion(:).numD])
        fig = [];
        figure('Units','centimeters','position',[20 20 16.97 24],'paperpositionmode','auto');
        for i=1:rows
            fig = gcf;
            figure(fig)
            subplot(rows,1,i)
            region = [pstart+(i-1)*(plength/rows) pstart+(plength/rows)*i];
            % plot the coverage
            plot(in(s).pos,movmean(in(s).cov,deletion(s).sldw),'Color',[0.8 0.8 0.8])
            hold on
            %plot the deletions
            maskde = (deletion(s).start >= region(1) & deletion(s).edge <= region(2));
            masknumde = sum(maskde);
            scatter(deletion(s).start(maskde),deletion(s).thresh*ones(masknumde,1),'d','filled','MarkerEdgeColor','m',...
                'MarkerFaceColor','m')
            scatter(deletion(s).edge(maskde),deletion(s).thresh*ones(masknumde,1),'d','filled','MarkerEdgeColor','b',...
                'MarkerFaceColor','b')
            %plot the duplications
            mask = (duplication(s).start >= region(1) & duplication(s).edge <= region(2));
            masknum = sum(mask);
            scatter(duplication(s).start(mask),duplication(s).thresh*ones(masknum,1),'d','filled','MarkerEdgeColor','m',...
                'MarkerFaceColor','m')
            scatter(duplication(s).edge(mask),duplication(s).thresh*ones(masknum,1),'d','filled','MarkerEdgeColor','b',...
                'MarkerFaceColor','b')
            %plot the thresholds and mean
            plot([region(1) region(2)],[deletion(s).mu deletion(s).mu],'k')
            hold on
            plot([region(1) region(2)],[deletion(s).thresh deletion(s).thresh],'k--')
            plot([region(1) region(2)],[duplication(s).thresh duplication(s).thresh],'k--')
            xlim([pstart+(i-1)*(plength/rows) pstart+(plength/rows)*i])
            ylim([0 1500])
            clear mask masknum maskde masksumde
        end
        subplot(rows,1,1)
        title([covsample(s) 'mapped to 168'])
    end
end
