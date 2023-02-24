function [gene,Gene_total,glib] = genes2fcn(geneMatrix,geneNames,nm_frmt,categories,cat_tmpl,varargin)
%Sorts a list of genes (gene names or gene IDs) into their respective
%categories. 

%INPUT
%genes -- list of  genes in matrix form, either as gene names OR as gene ID
%(BSU numbers)
%
%nm_frmt -- 1, gene ID (BSU number) 
%           2, gene name 
%
%categories -- list of all possible genes and categories (see SubtiWiki,
%gene categories as a reference)
%
%cat_tmpl -- the template used to interpret the gene categories file. Allows
%future users to use the script with different "categories" files
%
%OPTIONAL
%catlevel -- level/rank on which the categories should be analyzed--i.e,
%domain = 1, kingdom = 2, phylum =3, etc. This is an optional argument for
%the 'subtiwiki' template (cat_tmpl). If the argument is omitted, rank 1 is
%assumed.


paren = @(x, varargin) x(varargin{:});

%subtiwiki categories template
%find sheet name

if lower(cat_tmpl)=='subtiwiki'
    for k = 1:numel(varargin)
        switch lower(varargin{k})
            case 'sheet'
                Sheet = string(varargin{k+1});
        end
    end
    A = exist('Sheet','var');
    if A==0
        Sheet = '1';
        disp(['subtiwiki option requires the additional variable "sheet". ',...
            'followed by the sheet name. Assuming sheet name "1"']);
    end
    [~,~,genetable] = xlsread(categories,Sheet);
    cat_level = 1;
    nm_frmt = nm_frmt+1; %column in geneCategories spreadsheet
end

%Optional arguments
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'catlevel'
        cat_level = varargin{2};
        case 'sheet'
            if lower(cat_tmpl)=='subtiwiki'
                Sheet = varargin{2};
            end
        otherwise
        disp(['Unexpected option "', varargin{1},'" was ignored.'])
    end
    varargin(1:2)=[];
end

%%Gene type plots for CNPs
N = []; Ninfo=[];
for count = 1:size(geneMatrix,1)
    clear Info EYEnfo EYEnfox
    if isa(geneMatrix,'cell')
        Info = geneNames{count}';
    else
        Info = (geneNames(geneMatrix(count,:)>0))'; 
    end
    count2=0;u = 0;
    for i=1:size(Info,1)
        count2=count2+1;
        idx_gc = find(strcmp(Info(count2,1),genetable(:,nm_frmt))~=0);
        if isempty(idx_gc)
            EYEnfo(count2) = cellstr('Unknown'); %genes without functions (SubtiWiki list)
            %are marked as UNKNOWN
        elseif numel(idx_gc)>1
            EYEnfo(count2)=genetable(idx_gc(1),cat_level+3);
            for uu = 2:numel(idx_gc) %For genes with multiple functions, additional
                %(2nd, 3rd, etc.) functions are also added
                %to the EYEnfo list (at the end)
                EYEnfox(u+uu-1) = genetable(idx_gc(uu),cat_level+3);
            end
            u = u + uu -1;
        else
            EYEnfo(count2) = genetable(idx_gc,cat_level+3);
        end
    end
    if exist('EYEnfox','var')
        EYEnfo = cat(2,EYEnfo,EYEnfox);
    end
    [C{count},~,ic]=unique(Info(:,1));
    Cfreq{count} = histcounts(ic);
    if exist('EYEnfo','var')
        [G{count},~,ig]=unique(EYEnfo(:));
        Gfreq{count} = histcounts(ig);
    else
        G{count}=0;Gfreq{count}=0;
    end
    [~,~,Nin]=unique(Info(:,1));
    N=cat(1,N,cat(2,count*ones(length(Nin),1), Nin));
    Ninfo = cat(1,Ninfo,Info(:,1));
end

[glib, ~,iglib] = unique(genetable(:,cat_level+3));% 'Unknown'];
glib_cts = histcounts(iglib);
gene = zeros(numel(G),numel(glib));
for p = 1:numel(G)
    [~,idx_typ]=ismember(glib,G{p});
    for t = 1:numel(glib)
        if idx_typ(t) ~=0
            gene(p,t) = Gfreq{p}(idx_typ(t));% + gene(p,t);
        end
    end
end
glib_perc = glib_cts/4421;
Gene_total = gene./(cellfun(@sum,Cfreq(:))-cellfun(@(x) x(1,end),Gfreq(:))); %- ones(size(geneMatrix,1),1)*glib_perc; %gene./glib_cts;

%%%vvv%%%normalized based on how many genes were detected experimentally%%%vvv%%%
% gdem = ones(1,size(gene,2)).*cellfun(@sum,Cfreq)';
% Gene_total = gene./gdem;
%%%^^^%%%%%%^^^%%%%%%^^^%%%%%%^^^%%%%%%^^^%%%%%%^^^%%%%%%^^^%%%
%%%vvv%%%Only needed for frequency normalized per base pair%%%vvv%%%
% [~,BSUidx]=ismember(genetable(:,1),Nnames);
% for a = 1:max(iglib)
%     if sum(BSUidx(iglib==a))>0
%         temp = BSUidx(find(iglib==a & BSUidx~=0));
%         %the sum of the length of all genes in a given iglib category
%         glib_bp(a) = sum(gene168.L(Nnames_idx(temp) )); %assumes all genes are completely replaced
%     else
%         %some categories don't have any BSU matches
%         glib_bp(a)=0;
%     end
% end
% Gene_total = gene./glib_bp;
%%%^^^%%%%%%^^^%%%%%%^^^%%%%%%^^^%%%%%%^^^%%%%%%^^^%%%%%%^^^%%%

%Mean value over all samples in subset
mean(Gene_total,1)

figure(11)
bar(1:length(glib),paren(mean(Gene_total,1),1:numel(glib)));
hold on;
%calculate standard error of the mean
% if size(Gene_total,1)~=1
%     errL = paren(std(Gene_total),1:numel(glib));%/sqrt(numel(G));
%     %inlcude error measured from simulation
%     load GeneClass2Err.mat
%     errL = sqrt(errL.^2 + err.^2);
%     
%     %only plot error bars on non-zero mean values
%     errL(errL==0)=nan;
%     errorbar(paren(1:length(glib),~isnan(errL)),paren(mean(Gene_total,1),~isnan(errL)),...
%     paren(errL,~isnan(errL)),'ro','LineWidth',1.5);
% else
% %     load GeneClass2Err.mat
% %     errL = err;
% end

set(gca,'xtick',1:length(glib),'xticklabel',glib)
set(gca,'xticklabelrotation',65);
grid on
xlim([0 length(glib)+1])
ylabel('\Delta gene class frequency');
set(gca,'FontSize',12)
colormap(lines(1));
