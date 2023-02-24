function [deletion, duplication] = Cov2DelDub_fcn(covsample, covpath, covsuffix, artefacts, multimapper, recipsize, c, s,minL,delp,dupp,mincdel,mincdup)
%
%%%%%%%%%%%%%%%%   Cov2DelDup_fcn created by Isabel  %%%%%%%%%%%%%%%%%%%%%%%%
%
%
% What this script does: It finds all the deletions and duplications for
% given coverage data
%
%    input:
% -- is the coverage of the covsamples of interest; the artefacts to be excluded
% -- -- covsample, covpath, covsuffix as strings in cells
% -- -- if you want you can exclude artefacts and multimapper, but these
%       variables can also stay empty []
% -- and the parameters of choice
%    output:
% -- is a structure deletion e will pand duplication 


%% give the paths/names and set the parameters
if ~isempty(artefacts)
arte.path = artefacts;
fid = fopen(arte.path);
imp = textscan(fid,'%f %f %s %s','headerLines', 1);
fclose(fid);
arte.InS = imp{1};
arte.InE = imp{2};
arte.type = imp{3};
end
if ~isempty(multimapper)
%load mm regions and exclude them
mm.path = multimapper;
fid = fopen(mm.path);
imp = textscan(fid,'%f %f %f');
fclose(fid);
mm.InS = imp{1};
mm.InE = imp{2};
end

clrmap = jet(numel(covsample));

% variables 
refchr=recipsize; 
cutoff = c; % this is the factor c: cutoff= mu-c*sigma
sldw = s; % sliding window that smooths the cov prior to analysis
minL = minL; % minimal length of detected segment
del_para = delp; % this is the absolute value that del has to hit
min_consdel = mincdel; % del_para has to be hit this number of consecutive times
dup_para = dupp; % this is the factor to the mean that duplication has to hit
min_consdup = mincdup; % this is how often dub_para has to be hit consecutively

yes=1; no=0;
exclude = input('Do you want to exclude artefacts and multi mapper from the detectable regions? (yes/no) ');

%% preparations: 
% load the data
for i=1:numel(covsample)
    fid = fopen([covpath{:} covsample{i} covsuffix{:}]);
    imp = textscan(fid,'%s %f %f');
    fclose(fid);
    in(i).cov = imp{3};
    in(i).covsample = covsample{i};
    in(i).mu = mean(in(i).cov);
    in(i).sigma = std( in(i).cov);
    in(i).length = numel(in(i).cov);
    in(i).pos = 1:numel(in(i).cov);
    
    %plot the histograms as an overview
    figure(1)
    h(i)=histogram(in(i).cov,'FaceColor',clrmap(i,:),'EdgeColor',clrmap(i,:),'FaceAlpha',0.2,'BinWidth',5)
    hold on
    xlim([0 1000]) 
end

% make an artefacts mask for del and dup
artedel_mask = ones(refchr,1); artedup_mask = ones(refchr,1);
if ~isempty(artefacts)
delidx = find(ismember({arte.type{:}},'del'));
dupidx = find(ismember({arte.type{:}},'dup'));
end
if ~isempty(multimapper)
edel = [arte.InS(delidx) arte.InE(delidx); mm.InS mm.InE];
edup = [arte.InS(dupidx) arte.InE(dupidx); mm.InS mm.InE];
end
if ~isempty(artefacts) || ~isempty(multimapper)
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
    if isempty(exclude) || exclude == 1
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
    delmask2 = (movmean(in(i).cov,sldw) <= del_para);  
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
    deletion(m).covsample = covsample{i};
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
    if  isempty(exclude) || exclude == 1
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
    duplication(m).covsample = covsample{i};
    duplication(m).start = start_minL(maskmin0);
    duplication(m).edge = edge_minL(maskmin0);    
    % save other information
    duplication(m).numD = size(duplication(m).start,2);% - mockdel.numD;
    duplication(m).sumL = sum(duplication(m).edge - duplication(m).start +1);% - mockdel.refLsum;
end
%clean up if wanted
clear an dupmask dupmask2 edge edge_minL fid imp j m maskL maskmin0 ORI p stapel start start_minL sta sta2 start ed ed2 edge
