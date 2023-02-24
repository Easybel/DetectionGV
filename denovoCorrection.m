function [denovoC,varargout] = denovoCorrection(denovo,varannfolder,varannfiles,varargin)
%Correct the list of denovo mutations for known alignment errors. Also
%allow for the denovo matrix to include indels and SPIs.
%
%%%%%%%%%
%INPUTS
%denovo -- X by 1 cell of denovo mutation positions, for X samples.
%
%varannfolder -- folder where the annotated variant calls can be found (vcf
%files).
%
%varannfiles -- names of the individual annotated variant call files (vcf).
%The expected format is {'name1','name2',...}, and the default suffix is
%'_snpEff.vcf'
%
%%%%%%%%%%
%VARARGIN
%spi -- expects an X by 1 cell of SPI positions, where X is the number of
%samples
%
%filter1 -- filters out known INDEL alignment errors between Bsu168 and BsuW23.
%Requires a link to W23_indels.txt (tab separated POSITION, REF, ALT).
%
%filter2 -- filters out known artefacts between the recipient fasta
%reference and the experimental ancestral strain. Requires a link to
%artefacts_list_sort.txt (POSITION, REF, ALT, X of Y samples, GENENAME
%(when known)).
%
%region -- determines if the denovo mutation is inside or outside of a
%cluster. Requires the list of clusters, C.
%
%subset -- select a subset of the given denovo matrix. Expected format is
%[sample1 sample2 ... sampleN]. Matlab matrix operations are allowed--i.e.,
%1:10 and 1:2:10. VARANNFILES must have the same number of entries as
%subset.
%
%filesuffix -- suffix for the individual files. The default is '_snpEff.vcf'
%
%%%%%%%%%%%
%OUTPUTS
%denovoC -- corrected denovo mutation list. Format is an X by 1 cell, where 
%for each of the X samples, there is a list of mutation positions
%
%%%%%%%%%%%
%VARARGOUT
%R -- X by 1 cell containing binary information if the denovo mutation was
%inside (1) or outside (0) a cluster


%%%Load variables%%%

subset = 1:size(denovo,1);
filesuffix='_snpEff.vcf';
reg = 0; excpt1.N=0; excpt2.N=0; %Default, mutations are not grouped into regions (within/outside of clusters)
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
        case 'spi'
            SPI = varargin{2};
            for a = 1:numel(subset)
                temp{a} = sort(cat(2,denovo{subset(a)},SPI{subset(a)}));
            end
            clear denovo
            denovo = temp;
            clear temp
        case 'filesuffix'
            filesuffix=varargin{2};
        case 'subset'
            %previously addressed
        case 'region'
            reg = 1;
            Cluster = varargin{2};
            R = cell(numel(subset),1);
        case 'filter1'
            fid = fopen(varargin{2});
            EXCPT = textscan(fid,'%f %s %s','delimiter',' ');
            fclose(fid);
            excpt1.N=EXCPT{1}';
            excpt1.ref=EXCPT{2}';
            excpt1.alt=EXCPT{3}';
            clear EXCPT
            e1=1;
        case 'filter2'
            fid = fopen(varargin{2});
            EXCPT = textscan(fid,'%f %s %s %s %s','delimiter','\t');
            fclose(fid);
            excpt2.N=EXCPT{1}';
            excpt2.ref=EXCPT{2}';
            excpt2.alt=EXCPT{3}';
            clear EXCPT
            e2=1;
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
    filename = num2str(n);
    fid = fopen([varannfolder,varannfiles{subset(n)},filesuffix]);
    %Expected snpEff file format:
    %CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  20
    sample = textscan(fid,'%s %f %s %s %s %f %s %s %s %s','commentstyle','#');
    fclose(fid);
    clear info Info EYEnfo EYEnfox Jnfo
    info = sample{8}; place = sample{2}; ref = sample{4}; alt = sample{5};
    clear sample
    u = 0; count2=0;
    
    %Canidates for the corrected denovo matrix:
    A =find(ismember(place,denovo{n}) | strlength(ref)>1 | strlength(alt)>1);
    for i=1:numel(A)
        
    %Check if point mutations from the denovo list (opt., +SPI
    %list) or indels from the snpEff file need to be filtered using
    %artifact lists. If no artifact lists are given, filtering is
    %bypassed
            [memb, memb_idx]=ismember(place(A(i),1),excpt1.N);
            [memb2, memb2_idx]=ismember(place(A(i),1),excpt2.N);

            if memb_idx==0 && memb2_idx==0 %default; the mutation is filter approved
                FiltAppr = 1;
            elseif memb2_idx~=0 %excpt2 (requires alt base pair(s) to be identical)
                FiltAppr=~strcmp(excpt2.alt(memb2_idx),alt(A(i),1));
                if FiltAppr==0
                     countS(5,n)=countS(5,n)+1;
                end
            elseif memb_idx~=0 %&& FiltAppr==1 %excpt1 (requires alt base pair(s) to be identical)
                                              %If there are multiple alt,
                                              %the first mutation is used
                if contains(excpt1.alt(memb_idx),',') || contains(alt(A(i),1),',')
                    clear splt_alt splt_indelalt
                    splt_alt = strsplit(char(alt(A(i),1)),',');
                    splt_indelalt = strsplit(char(excpt1.alt(memb_idx)),',');
                    if sum(ismember(splt_alt,splt_indelalt))>0
                        FiltAppr = 0;
                    end
                else
                    FiltAppr=~strcmp(excpt1.alt(memb_idx),alt(A(i),1));
                end
            end
        
            if FiltAppr==1
                count2=count2+1;
                if reg ==1 %record in cluster/outside of cluster information if option 'region' is given
                    if sum(place(A(i),1)>=min(Cluster{n}(1:2,:)) & place(A(i),1)<=max(Cluster{n}(1:2,:))) ==1
                        R{n,1}(1,count2) = 1;
                    else 
                        R{n,1}(1,count2) = 0;
                    end
                end
                clear temp
                denovoC{n,1}(1,count2)=place(A(i),1);
            end
    end
    denovoC{n,1} = sortrows(denovoC{n,1});
end


if nargout==2
    varargout{1} = R;
end






