function [Cdist,Adist,Ident, varargout] = lengthNident(recipgenome,recipsize,masterlist,Cluster,varargin)
%Calculates the maximum, minimum, and mean CNP sizes (from Cluster{m}) along
%with the identies of each CNP.
%
%%%%%%%%%
%INPUTS
%recipgenome - recipient genome auxiliary genes. Five tab separated
%columns. (1) generic gene reference (text), (2) start position, (3) end
%position, (4) length (aka (3)-(2) plus one), (5) gene description
%---'/scratch/bobbypickett/DirectsW23/168genes.csv'
%
%recipsize -- recipient genome size (bp)
%
%masterlist -- master list for the recipient/donor pair. Two tab separated
%columns. (1) position (2) alternate nucleotide
%
%Cluster -- X by 1 cell of cluster positions. For each of the X samples in
%cluster, there are 4 rows: start position (Row1), end position (Row2),
%length (Row3), number of missing master list SNPs (Row4)
%
%%%%%%%%%
%VARARGIN
%diagnostics -- (default OFF (0)) this option causes the function to produce plots related to
%determine how effective the clustering is (# missing SNPs, etc.).
%
%mdist -- measure the maximum length of a cluster from the first
%not-detected master list SNP before and after the cluster
%
%ident2 -- measure the identity of cluster segments taking into account
%missing SNPs
%
%OUTPUTS
%Cdist -- length of cluster elements from first detected master list SNP to 
%last detected master list SNP. An X by 1 cell containing lengths.
%
%Adist -- average length of cluster elements, averaging the detected length
%(Cdist) and the maximum length (Mdist). An X by 1 cell containing lengths.
%
%Ident -- identity of the cluster segments. An X by 1 cell (decimals, not
%percentages)
%
%%%%%%%%%%%
%VARARGOUT
%from 'diagnostics'
%Fdist -- distance from cluster to next position on the master list
%nCNP -- structure with stat (measure using cluster sliding window, see
%code below) and dist (distance between clusters)
%ib -- structure with numb (of In Between cluster master list positions),
%pos (their positions), dist (the distance between each of them), and ident
%(identity if all master list SNPs in the region were present)

%from 'mdist'
%Mdist -- maximum length of a cluster from the first not-detected 
%master list SNP before and after the cluster. X by 1 cell with lengths

%from 'ident2'
%Ident2 -- identity of cluster segments taking into account missing SNPs. X
%by 1 cell (fractions not percentages)


%%%Load variables and optional arguments%%%%%%%%%
diag = 0;
%Optional inputs
k=0;
while k<numel(varargin)
    k=k+1;
    switch lower(varargin{k})
        case 'diagnostics'
            diag = 1; template = cell(1:numel(size(Cluster,2),1));
            nCNP = struct('stat',[],'dist',[]);
            ib = struct('numb',[],'pos',[],'dist',[],'ident',[]);
            %             ib_miss.pos = cell(1:numel(size(Cluster,2),1));
            %             ib_miss.dist = cell(1:numel(size(Cluster,2),1));
            acpwin = varargin{k+1};
            if isa(varargin{k+1},'double')
                k = k+1;
            else
                disp(['Unexpected option after ''diagnostics'': ', varargin{k+1},...
                    '. Expected entry was CLUSTER WINDOW SIZE. Dignostics option ignored.'])
                diag = 0;
            end
        case 'mdist'
        case 'ident2'
        otherwise
            disp(['Unexpected option "', varargin{k},'" was ignored.'])
    end
end

%Load recipient/donor specific master list
fid = fopen(masterlist);
imp = textscan(fid,'%f %s');
fclose(fid);
refchr.pos=imp{1};
refchr.atcg=imp{2};
clear imp

%Define variables
template = cell(1:numel(size(Cluster,2),1));
Cdist=template; Mdist=template; Adist=template; Fdist=template;
Ident=template; Ident2=template;
% LCNP=template;
% N = zeros(numel(size(Cluster,2)),ceil(recipsize/1e5)); TF=nan(numel(size(Cluster,2)),4);


for m = 1:size(Cluster,1)
    if ~isempty(Cluster{m})
        %Cluster size from first to last detected SNP
        Cdist{m}=max(Cluster{m}(1:2,:)) - min(Cluster{m}(1:2,:)) +1; %+1 to include starting SNP
        %Identity (for integrated segment (CNP) compared to recipient genome
        Ident{m} = (Cdist{m} - Cluster{m}(3,:))./Cdist{m};
        %Ident2 is the identity assuming the entire donor strand was
        %recombined [aka adding back the missing SNPs (Cluster{}(4,:))]
        Ident2{m} = (Cdist{m} - sum(Cluster{m}(3:4,:)) )./Cdist{m};
        
        %Measure maximum import lengths from minimum lengths (Cdist)
        %         MCluster{m} = zeros(size(Cluster{m})); %Fdist=zeros(2,size(Cluster{m},2));
        for p = 2:size(Cluster{m}(1,:),2)-1
            %first possible masterlist SNP, right of a CNP minus first
            %possible masterlist SNP, left of a CNP
            hitL= refchr.pos(find(refchr.pos < min(Cluster{m}(1:2,p)),1,'last'));
            hitR= refchr.pos(find(refchr.pos > max(Cluster{m}(1:2,p)),1,'first'));
            Mdist{m}(1,p)=hitR - hitL + 1; %Count the position one is on
            Fdist{m}.L(p) = min(Cluster{m}(1:2,p)) - hitL; %Fdist is the distance to
            %the next master list SNP. The matrix is only used with the
            %'diagnostics' option
            Fdist{m}.R(p) = hitR - max(Cluster{m}(1:2,p));
        end
        %First and last clusters could jump over origin
        %First cluster
        hitL = refchr.pos(find(refchr.pos < min(Cluster{m}(1:2,1)),1,'last'));
        hitR= refchr.pos(find(refchr.pos > max(Cluster{m}(1:2,1)),1,'first'));
        if isempty(hitL)
            Mdist{m}(1,1)= recipsize - refchr.pos(end) +...
                refchr.pos(find(refchr.pos > max(Cluster{m}(1:2,1)),1,'first')) + 1;
            Fdist{m}.L(1) = recipsize - refchr.pos(end) + min(Cluster{m}(1:2,1));
            Fdist{m}.R(1) = hitR - max(Cluster{m}(1:2,1));
        else
            Mdist{m}(1,1) = hitR - hitL + 1;
            Fdist{m}.L(1) = min(Cluster{m}(1:2,1)) - hitL;
            Fdist{m}.R(1) = hitR - max(Cluster{m}(1:2,1));
        end
        %Last cluster
        hitR = refchr.pos(find(refchr.pos > max(Cluster{m}(1:2,end)),1,'first'));
        hitL = refchr.pos(find(refchr.pos < min(Cluster{m}(1:2,end)),1,'last'));
        if Cluster{m}(2,end)<1 || Cluster{m}(2,end)>recipsize
            hitR = refchr.pos(find(refchr.pos > min(mod(Cluster{m}(1:2,end),recipsize)),1,'first'));
            Mdist{m}(1,end+1) = recipsize + hitR - hitL + 1;
            Fdist{m}.L(end+1) = mod(min(Cluster{m}(1:2,end)),recipsize) - hitL;
            Fdist{m}.R(end+1) = hitR - min(mod(Cluster{m}(1:2,end),recipsize));
        elseif isempty(hitR)
            Mdist{m}(1,end+1) = recipsize - hitL + refchr.pos(1) + 1;
            Fdist{m}.L(end+1) = mod(min(Cluster{m}(1:2,end)),recipsize) - hitL;
            Fdist{m}.R(end+1) = recipsize + refchr.pos(1) - max(Cluster{m}(1:2,end));
        else
            Mdist{m}(1,end+1) = hitR - hitL +1;
            Fdist{m}.L(end+1) = min(Cluster{m}(1:2,end)) - hitL;
            Fdist{m}.R(end+1) = hitR - max(Cluster{m}(1:2,end));
        end
        
        %"Average" distance -- average of Mdist and Cdist
        Adist{m} = mean([Cdist{m}; Mdist{m}]);
        
        
        
        %Diagnostic calculations/plots if the option ('diagnostics',1) is
        %selected
        if diag == 1
            %nCNP = non-CNP
            %nCNP(1) is the non Cluster{m} region immediately after CNP(1)
            
            %If next nearest SNP was <acpwin bp away from the Cluster{m}
            %aka =<(acpwin-2), because the last and first SNPs count as one of the 200 bp
            nCNP(m).stat(Fdist{m}.L(2:end)>=acpwin-1 & Fdist{m}.R(1:end-1)>=acpwin-1) = 0;%0= >200bp
            nCNP(m).stat(Fdist{m}.L(2:end)<acpwin-1 & Fdist{m}.R(1:end-1)>=acpwin-1) = 1;%1= <200bp one side
            nCNP(m).stat(Fdist{m}.L(2:end)>=acpwin-1 & Fdist{m}.R(1:end-1)<acpwin-1) = 1;
            nCNP(m).stat(Fdist{m}.L(2:end)<acpwin-1 & Fdist{m}.R(1:end-1)<acpwin-1) = 2;%2= <200bp both sides
            nCNP(m).stat(2:end)=nCNP(m).stat(1:end-1); %matlab indexs improperly above
            nCNP(m).stat(Fdist{m}.L(1)>=acpwin-1 & Fdist{m}.R(end)>=acpwin-1) = 0;
            nCNP(m).stat(Fdist{m}.L(1)<acpwin-1 & Fdist{m}.R(end)>=acpwin-1) = 1;
            nCNP(m).stat(Fdist{m}.L(1)>=acpwin-1 & Fdist{m}.R(end)<acpwin-1) = 1;
            nCNP(m).stat(Fdist{m}.L(1)<acpwin-1 & Fdist{m}.R(end)<acpwin-1) = 2;
            
            %Distance BETWEEN CNPs
            nCNP(m).dist(1:size(Cdist{m},2)-1) = Cluster{m}(1,2:end)-Cluster{m}(2,1:end-1) -1;
            nCNP(m).dist(size(Cdist{m},2))   = recipsize - Cluster{m}(2,end) + Cluster{m}(1,1) -1;
                        
            for w = 1:size(Cdist{m},2)-1
                %ib -- in between CNPs, number of potential SNPs,
                %positions, distances, and identity
                ib(m).numb(w) = numel(find(refchr.pos<Cluster{m}(1,w+1) & refchr.pos>Cluster{m}(2,w)));
                ib(m).pos = cat(1,ib(m).pos,refchr.pos(refchr.pos<Cluster{m}(1,w+1) & refchr.pos>Cluster{m}(2,w)));
                ib(m).dist = cat(2,ib(m).dist,diff(refchr.pos(refchr.pos<Cluster{m}(1,w+1) &...
                    refchr.pos>Cluster{m}(2,w)))');
                
                %             if diag == 1
                %                 %Pull up matches (MATCH.m) to look for the distance between
                %                 %missing SNPs
                %                 match = MATCH(sample,refchr);
                %
                %                 %Distance between missing SNPs
                %                 ib_miss(m).pos = cat(1,ib_miss(m).pos,setdiff(refchr.pos(refchr.pos>=Cluster{m}(1,w) &...
                %                     refchr.pos<=Cluster{m}(2,w)),match));
                %                 ib_miss(m).dist = cat(2,ib_miss(m).dist,diff(setdiff(refchr_SNPs(refchr_SNPs>=Cluster{m}(1,w) &...
                %                     refchr_SNPs<=Cluster{m}(2,w)),match)'));
                %                 if ib(w) ==0
                %                     nCNP_stat(w)=-1;% -1=no SNPs in nCNP
                %                 end
                %             end
                
                %Identity of BETWEEN CNP segements
                ib(m).ident(w) = ( nCNP(m).dist(w) - ib(m).numb(w) )/  nCNP(m).dist(w);
                
            end
            %         ib(numel(cdist)) = numel(find(refchr_SNPs<Cluster{m}(1,1) |...
            %             refchr_SNPs>Cluster{m}(2,numel(cdist))));
            %         ib_pos_d = cat(2,ib_pos_d,diff(refchr_SNPs(refchr_SNPs<Cluster{m}(1,1) |...
            %             refchr_SNPs>Cluster{m}(2,numel(cdist))))');
            ib(m).ident(size(Cdist{m},2)) = ( nCNP(m).dist(size(Cdist{m},2))-ib(m).numb(end))/...
                nCNP(m).dist(size(Cdist{m},2));
            if Cluster{m}(1,size(Cdist{m},2))>Cluster{m}(2,size(Cdist{m},2))
                %             LCNPt(m_idx,numel(cdist)) = numel(find(refchr_SNPs>=Cluster{m}(1,numel(cdist))))+...
                %                 numel(find(refchr_SNPs<=Cluster{m}(2,numel(cdist))));
                ib(m).pos = cat(1,ib(m).pos,refchr.pos(refchr.pos<Cluster{m}(1,1) &...
                    refchr.pos>Cluster{m}(2,size(Cdist{m},2))));
            else
                %             LCNPt(m_idx,numel(cdist)) = numel(find(refchr_SNPs>=Cluster{m}(1,numel(cdist)) &...
                %                 refchr_SNPs<=Cluster{m}(2,numel(cdist))));
                ib(m).pos = cat(1,ib(m).pos,refchr.pos(refchr.pos<Cluster{m}(1,1) |...
                    refchr.pos>Cluster{m}(2,size(Cdist{m},2))));
            end
            %         LCNP(m_idx,numel(cdist)) = LCNPt(m_idx,numel(cdist)) - Cluster{m}(3,numel(cdist));
            %         ib_miss(m).pos = cat(1,ib_miss(m).pos,setdiff(refchr_SNPs(refchr_SNPs>=Cluster{m}(1,numel(cdist)) &...
            %             refchr_SNPs<=Cluster{m}(2,numel(cdist))),match));
            %         ib_miss(m).dist = cat(2,ib_miss(m).dist,diff(setdiff(refchr_SNPs(refchr_SNPs>=Cluster{m}(1,numel(cdist)) &...
            %             refchr_SNPs<=Cluster{m}(2,numel(cdist))),match)'));
            %         %Litmus test for nCNPs (# SNPs, length, =>< 200bp)
            %         LnCNP{m_idx,:}=[ib; FCNPdist(m_idx,1:numel(cdist)); nCNP_stat];
            %         TF(m_idx,:) = [sum(Cluster{m}(3,:)), nansum(LCNP(m_idx,:)),...
            %             numel(ib_pos)-numel(intersect(ib_pos,match)),numel(intersect(ib_pos,match))];
            %         SPI{m_idx,:}=[intersect(ib_pos,match); setdiff(Sample,refchr_SNPs)];
            %CNP Litmus
            
            %             Fsumdist(m_idx,1:numel(cdist))= sum(Fdist);
            %             figure(4)
            % %             bins = 1:5:501; %Left edge is included in the bin
            %             histogram(FLdist,bins,'normalization','probability'); hold on;
            %             histogram(FUdist,bins,'normalization','probability')
            %             xlim([0 501]);
            %
            %             figure(1)
            %             edges = 10.^(0:0.1:4);
            %             histogram(adist,edges)%,'normalization','probability')
            %             set(gca,'XScale','log'); grid on; hold on;
        end
    end
end

%Transpose Cdist, Adist, Ident
Cdist = Cdist';
Adist = Adist';
Ident = Ident';

k=0;
while ~isempty(varargin)
    
    switch lower(varargin{1})
        case 'diagnostics'
            if diag==1
                varargin(2)=[];
                varargout{k+1}=Fdist';
                varargout{k+2}=nCNP;
                varargout{k+3}=ib;
                k = k+3;
            end
        case 'mdist'
            varargout{k+1} = Mdist';
            k = k+1;
        case 'ident2'
            varargout{k+1} = Ident2';
            k = k+1;
    end
    varargin(1)=[];
end

