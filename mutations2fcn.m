function [Xno,Gene_total,Jean] = mutations2fcn(denovoC,varannfolder,varannfiles,varargin)
%Sort de novo mutations based upon their predicted impact
%
%INPUTS
%region --
%
%OUTPUTS
%
%OPTIONAL
%

subset = 1:size(denovoC,1);
filesuffix='_snpEff.vcf';
reg = 0; fcn=0; %Defaults
cat_level = 1;
paren = @(x, varargin) x(varargin{:});

%Optional arguments
k=1;
while k<numel(varargin)
    switch lower(varargin{k})
        case 'subset'
            if isa(varargin{k+1},'double') && max(varargin{k+1}) <= subset(end) &&...
                    numel(varargin{k+1})<=numel(varannfiles)
                subset = varargin{k+1};
                k = k+2;
            else
                disp(['Unexpected option after ''subset'': ', varargin{k+1},...
                    '. Expected entry is a double. Subset option ignored.'])
            end
        otherwise
            k = k+2;
    end
end
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'filesuffix'
            filesuffix=varargin{2};
        case 'subset'
            %previously addressed
        case 'region'
            reg = 1;
            R = varargin{2};
        case 'fcn'
            fcn = 1;
            [~,~,genetable] = xlsread(varargin{2},'geneCategories (1)');
        otherwise
            disp(['Unexpected option "', varargin{1},'" and subsequent option was ignored.'])
    end
    varargin(1:2)=[];
end


N = []; Ninfo=[];
countS = zeros(5,numel(subset)); %row1=upstream inserts %row2=upstream deletions
%row3=upstream indels matching W23 %row4=INTRAgenic
%indels matching W23

for n = 1:numel(subset)
    clear sample
    fid = fopen([varannfolder,varannfiles{subset(n)},filesuffix]);
    %Expected snpEff file format:
    %CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  20
    sample = textscan(fid,'%s %f %s %s %s %f %s %s %s %s','commentstyle','#');
    fclose(fid);
    clear info Info EYEnfo EYEnfox Jnfo
    info = sample{8}; place = sample{2}; ref = sample{4}; alt = sample{5};
    clear sample
    [~,A] = intersect(place,denovoC{n});
    u = 0;
    for i=1:numel(A)
        clear temp
        
        temp=strsplit(char(info(A(i),1)),'|');
        %diagnostics
        %                 if strlength(ref(A(i),1))>1  && temp(2)=="upstream_gene_variant"
        %                     countS(1,i) = countS(1,i)+1;
        %                     if ismember(place(A(i),1),indel_N) %see if indels match W23/168 mapping
        %                         countS(3,i)=countS(3,i)+1;
        %                     end
        %                 elseif strlength(alt(A(i),1))>1 && temp(2)=="upstream_gene_variant"
        %                     countS(2,i)=countS(2,i)+1;
        %                     if ismember(place(A(i),1),indel_N) %see if indels match W23/168 mapping
        %                     countS(3,i)=countS(3,i)+1;
        %                     end
        %                 elseif (strlength(ref(A(i),1))>1 || strlength(alt(A(i),1))>1) &...
        %                         (temp(2)=="frameshift_variant" || temp(2)=="disruptive_inframe_deletion" ||...
        %                         temp(2)=="disruptive_inframe_insertion" || temp(2)=="inframe_insertion" ||...
        %                         temp(2)=="inframe_deletion")
        %                     if ismember(place(A(i),1),indel_N) %see if indels match W23/168 mapping
        %                     countS(4,i)=countS(4,i)+1;
        %                     end
        %                 end
        if temp(3)== "ERROR_OUT_OF_CHROMOSOME_RANGE"
            Info(i,1)={'dnaA'};%Indels at the end of the chromosome are
            %upstream of the first gene, dnaA
            Info(i,2)={'frameshift variant'};
            Info(i,3)={'BSU00010'};
            %%%%%%%%%%%%%vvvvvvvvvvvvvvvvvvvvvvvvvvvv%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif contains(temp(2),'&')
            if contains(temp(2),'start') ||contains(temp(2),'initiator')
                Info(i,2)={'start'};
            elseif contains(temp(2),'stop')
                Info(i,2)={'stop'};
            elseif contains(temp(2),'frameshift_variant')
                Info(i,2)={'frameshift variant'};
            elseif contains(temp(2),'disruptive_inframe_insert')
                Info(i,2)={'disruptive_inframe_insert'};
            else
                Info(i,2)=paren(split(temp(2),'&'),1);
            end
            Info(i,1)=strrep(temp(4),'_',' ');
            Info(i,3)=temp(5);
            %%%%%%%%%%%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            Info(i,1)=strrep(temp(4),'_',' ');
            Info(i,2)=strrep(temp(2),'_',' ');
            Info(i,3)=temp(5);
        end
        idx_gc = find(strcmp(Info(i,3),genetable(:,2))~=0);
        if isempty(idx_gc)
            EYEnfo(i) = cellstr('Unknown'); %genes without functions (SubtiWiki list)
            %are marked as UNKNOWN
            Jnfo(i) = cellstr('Non-essential genes');
        elseif numel(idx_gc)>1
            EYEnfo(i)=genetable(idx_gc(1),4);
            for uu = 2:numel(idx_gc) %For genes with multiple functions, additional
                %(2nd, 3rd, etc.) functions are also added
                %to the EYEnfo list (at the end)
                EYEnfox(u+uu-1) = genetable(idx_gc(uu),4);
            end
            u = u + uu -1;
            
        else
            EYEnfo(i) = genetable(idx_gc,4);
        end
        if contains('Essential genes',genetable(idx_gc(:),5))
            Jnfo(i) = cellstr('Essential genes');
        else
            Jnfo(i) = cellstr('Non-essential genes');
        end
    end
    %         end
    
    if exist('EYEnfox')
        EYEnfo = cat(2,EYEnfo,EYEnfox);
    end
    if ~exist('Info','var')
        C{n} = {'filler'};
        Cfreq{n}=1;
        G{n} = {'filler'};
        Jfreq{count} = [0 0];
    else
        [C{n},~,ic]=unique(Info(:,2));
        Cfreq{n} = histcounts(ic);
        [G{n},~,ig]=unique(EYEnfo(:));
        Gfreq{n} = histcounts(ig);
        [J{n},~,ij]=unique(Jnfo(:),'stable');
        Jfreq{n} = [length(find(not (cellfun('isempty',(strfind(Jnfo,"Essential genes")))))) length(find( (cellfun('isempty',(strfind(Jnfo,"Essential genes"))))))]; %histcounts(ij);
        [Nnames,~,Nin]=unique(Info(:,1));
        N=cat(1,N,cat(2,n*ones(length(Nin),1), Nin));
        Ninfo = cat(1,Ninfo,Info(:,1));
    end
end

clear Nnames
[Nnames,~,in]=unique(Ninfo);
Tmatch = zeros(size(C,2),size(Nnames,1));
for a = 1:size(C,2)
    f=find(N(:,1)==a);
    Tmatch(a,in(f))=1;
end
Tmatch2 = array2table(Tmatch);
Nnames2=unqNames(Nnames);
SampleNames = cellstr(string(subset)');
Tmatch2.Properties.VariableNames = Nnames2;
Tmatch2.Properties.RowNames = SampleNames;
% Tmatch = cell2table(cat(2,SampleNames,cat(1,Nnames2,num2cell(Tmatch))));
%writetable(Tmatch2,'GenesOfDeNovoMutations_cy15.xlsx','writevariablenames',true,'writerownames',true)


% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %gene type count
if fcn ==1
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
    Gene_total = gene./glib_cts;
        
    %Mean value over all samples in subset
%     mean(Gene_total)
    
    figure(11)
    bar(1:length(glib),paren(mean(Gene_total),1:numel(glib)));
    hold on;
    %calculate standard error of the mean
    err = paren(std(Gene_total),1:numel(glib))/sqrt(numel(G));
    %only plot error bars on non-zero mean values
    err(err==0)=nan;
    errorbar(paren(1:length(glib),~isnan(err)),paren(mean(Gene_total),~isnan(err)),...
        paren(err,~isnan(err)),'ro','LineWidth',1.5);
    set(gca,'xtick',1:length(glib),'xticklabel',glib)
    set(gca,'xticklabelrotation',65);
    grid on
    xlim([0 length(glib)+1])
    ylabel('Gene class frequency');
    set(gca,'FontSize',12)
    colormap(lines(1));
    
end

%%
%%%%%%%%%%%%%vvvvvvvvvvvvvvvvvvvvvvvvvvvv%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SNP type count
lib = {'synonymous variant','missense variant', ...
    'stop retained variant', 'stop gained','stop lost','stop',...
    'initiator codon variant','start retained', 'start lost','start' ...
    'frameshift variant',...
    'disruptive inframe deletion', 'disruptive inframe insertion', 'inframe insertion',  'inframe deletion', ...
    'upstream gene variant',...
    'intragenic variant','downstream gene variant',... %"catch all" for misfit variants
    'protein protein contact','structural interaction variant','rare amino acid variant', 'intergenic variant','coding sequence variant'};
%last 5 entries are flaged, should not be found
%%%%%%%%%%%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xno = zeros(numel(C),numel(lib));
for p = 1:numel(C)
    
    [~,idx_typ]=ismember(lib,C{p});
    for t = 1:numel(lib)
        if idx_typ(t) ~=0
            xno(p,t) = Cfreq{p}(idx_typ(t)) + xno(p,t);
        end
    end
    
end

lib_abv = {'syn mut','nonsyn mut', ...
    'stop mut', 'start mut', ...
    'indel missense', 'indel inframe', ...
    'upstream mut','intragen mut'};
%%%%%%%%%%%%%vvvvvvvvvvvvvvvvvvvvvvvvvvvv%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xno = [xno(:,1) xno(:,2) sum(xno(:,3:6),2) sum(xno(:,7:10),2)...
    xno(:,11) sum(xno(:,12:15),2) xno(:,16) sum(xno(:,17:18),2)];
%INTERgenic mutations are not displayed
flag=0;
if sum(sum(xno(:,19:23),2))~=0
    flag = 1;
end
%%%%%%%%%%%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gdem = ones(1,size(gene,2)).*sum(gene,2);
% Gene = gene./1;%gdem;
dem = ones(1,size(Xno,2)).*sum(Xno,2);
Xno_norm = Xno./dem;

figure(2)
bar(1:size(Xno,1),Xno(:,1:numel(lib_abv)),'stacked')
%ylim([0 1.05]);
legend(lib_abv,'location','northeastoutside')
ylabel('No. of occurances');
%ylabel('Total number of mutations');
grid on
%set(gca,'xtick',1:14,'xticklabel',{'W1','W2','W4','W5','W6','W7','W8','cy21'})
set(gca,'xticklabelrotation',0);
clr = colormap([lines(7); colorcube(4)]);
colormap(clr([7:-1:1,10],:))
xlim([0 8]);
if flag ==1
    title('INTERgen mut./p-p interaction/rare AA','Color','r')
    set(gca,'Color','r')
end
set(gca,'FontSize',12)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Essential genes

%find essential gene category
jean = cell2mat(Jfreq');
jean_r = sqrt(jean);
jdem = [255 4421-255]; %No. of essential genes, no. of non-essential genes
Jean = jean./sum(jean,2);
Jean_r = jean_r./jdem;
% Jean_total = jean./[glib_cts(hit) sum(glib_cts([1:hit-1,hit+1:end]))];

figure(12)
% boxplot(Jean(:,2)');
% set(gca,'xtick',1:8,'xticklabel',{'W1','W2','W3','W4','W5','W6','W7','W8'})
% hold on
plot(9,nanmean(Jean(:,2)'),'o');
hold on
errorbar(9,nanmean(Jean(:,2)'),nanstd(Jean(:,2)'))
%set(gca,'xticklabelrotation',35);
% legend('W1','W2','W3','W4','W5','W6','W7','W8','location','northwest');
%'no DNA','wt DNA','B.moj DNA','location','best');
grid on
xlim([0 10]); ylim([0.85 1])
plot([0 10],[(4421-255)/4421 (4421-255)/4421],'k--')
%ylabel('No. of mutations');
ylabel('Frequency of mutated-gene class');
colormap((lines(1)));

%%
%Plot of mutated genes
% rows=5;
% %Plot A or B (it is split over two images)
% pl = 1; %A = 1, B=2
% f10=figure('Units','centimeters','position',[20 20 16.97 24],'paperpositionmode','auto');
% k = [1:8];
% for b = 1:rows
%     subplot(rows,1,b)
%     hold on;
%     imagesc(flipud(~isnan(Tmatch)))
%     caxis([0 1])
% end
% ell = floor(length(Tmatch)/(2*(rows - (1/8))));
% Nnames_str=cellstr(Nnames);
% if pl==1
%     for b = 1:rows
%         subplot(rows,1,b)
%         xlim([0+ell*(b-1)+.5,ell*b+.5]); grid on;
%         ylim([0.5 numel(k)+.5]); ax = gca;
%         ax.YTick=[1:numel(k)];
%         ax.YTickLabel={'W8','W7','W6','W5','W4','-- ','W2','W1'};
%         ax.XTick=[0+ell*(b-1)+1:ell*b];
%         ax.XTickLabel=Nnames_str(0+ell*(b-1)+1:ell*b);
%         set(ax,'xticklabelrotation',90,'FontSize',12);
%     end
%     subplot(rows,1,rows)
%     xlim([0+ell*(b-1)+0.5,.5+length(Tmatch)/2]);
%
% else %pl=2
%     for b = 1:rows
%         subplot(rows,1,b)
%         xlim([0+ell*(b-1)+.5+length(Tmatch)/2,ell*b+.5+length(Tmatch)/2]); grid on;
%         ylim([0.5 numel(k)+.5]); ax = gca;
%         ax.YTick=[1:numel(k)];
%         ax.YTickLabel={'W8','W7','W6','W5','W4','--','W2','W1'};
%         if ell*b+length(Tmatch)/2>length(Nnames)
%             ax.XTick=[0+ell*(b-1)+1+length(Tmatch)/2:length(Tmatch)];
%             ax.XTickLabel=Nnames_str(0+ell*(b-1)+1+length(Tmatch)/2:length(Tmatch));
%         else
%             ax.XTick=[0+ell*(b-1)+1+length(Tmatch)/2:length(Tmatch)];
%             ax.XTickLabel=Nnames_str(0+ell*(b-1)+1+length(Tmatch)/2:ell*b+length(Tmatch)/2);
%         end
%         set(ax,'xticklabelrotation',90,'FontSize',12);
%     end
%     subplot(rows,1,rows)
%     xlim([0+ell*(b-1)+0.5+length(Tmatch)/2,length(Tmatch)+.5]);
%
% end
% subplot(rows,1,rows)
% colorbar('eastoutside','TickLabels',{'Not affected','Affected'},'Ticks', [0.25 0.75]);
% %cmap = colormap(gray);
% cmap = [.9683*ones(1,3); 0 0 0];
% colormap(cmap);
% subplot(rows,1,3)
% ylabel('W23 replicates, cycle 21');
%
%


