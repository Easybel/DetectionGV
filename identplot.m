function [h1,varargout] = identplot(Ident,Cdist,varargin)
%Scatter plot of identities
%
%INPUT
%Ident -- X by 1 cell array with identity values for each sample (as
%fractions, not percentage)
%
%Adist -- X by 1 cell array of average CNP lengths
%
%VARARGIN
%subset -- allows one to plot only a selection of the Ident cell. The
%following formating is expected: [1,2, ... N]. Matlab operations are
%allowed--i.e., 1:10 and 1:2:10.
%
%xlimit -- fix x range on plot
%
%ylimit -- fix y range on plot
%
%ave -- genome wide average sequence identity. The de
%
%%%%%%%%
%OUTPUT
%h1 -- handle to identity plot
%
%VARARGOUT
%option 'ave':
%ident -- structure containing two matrices .p (<100 bp segments) and .m
%(>=100bp segments). These statistics matrices contain mean (1), std (2),
%numel (3), p-value (4).


paren = @(x, varargin) x(varargin{:});
subset = 1:size(Ident,1);

%Optional variables
k=0; xlimit = 0; ylimit = 0;
while k<numel(varargin)
    k=k+1;
    switch lower(varargin{k})
        case 'subset'
            if isa(varargin{k+1},'double') && max(varargin{k+1}) <= subset(end)
                subset = varargin{k+1};
                k = k+1;
            else
                disp(['Unexpected option after ''subset'': ', varargin{k+1},...
                    '. Expected entry is a double. Subset option ignored.'])
            end
        case 'xlimit'
            if isa(varargin{k+1},'double')
                xlimit = varargin{k+1};
                k = k+1;
            else
                disp(['Unexpected option after ''xlimit'': ', varargin{k+1},...
                    '. Expected entry is a double. xlimit option ignored.'])
            end
        case 'ylimit'
            if isa(varargin{k+1},'double')
                ylimit = varargin{k+1};
                k = k+1;
            else
                disp(['Unexpected option after ''ylimit'': ', varargin{k+1},...
                    '. Expected entry is a double. ylimit option ignored.'])
            end
        case 'ave'
            if isa(varargin{k+1},'double')
                identAve = varargin{k+1};
                k = k+1;
            else
                disp(['Unexpected option after ''ave'': ', varargin{k+1},...
                    '. Expected entry is a double. ave option ignored.'])
            end
        otherwise
            disp(['Unexpected option "', varargin{k},'" was ignored.'])
            
    end
end

clr = (parula(numel(subset)));



h=figure(1);
h1=histogram(paren(reshapeBOX(Cdist',subset),':',1));
if xlimit==0
    xlimit = h1.BinLimits;    
end
if ylimit==0
    h2=histogram(paren(reshapeBOX(Ident',subset),':',1));
    ylimit = h2.BinLimits;
end
close(h)
clear h1 h2

h1=figure(1);

subplot(4,4,[14:16])
histogram(paren(reshapeBOX(Cdist',subset),':',1),'FaceColor',lines(1),'BinLimits',xlimit); 
hold on; grid on;
xlabel('Segment length (bp)'); ylabel('Freq.')
colormap(jet)
xlim([xlimit])

subplot(4,4,[1 5 9])
histogram(paren(reshapeBOX(Ident',subset),':',1),'FaceColor',lines(1),'BinLimits',ylimit );
camroll(90)
hold on; grid on;
set(gca,'XAxisLocation','top'); 
xlabel('Identity'); ylabel('No. events')
xlim([ylimit]);
A = exist('identAve');
if A~=0
    if identAve>1
        identAve = identAve/100;
    end
    histylimits=ylim;
    hold on; plot([identAve identAve],[histylimits],'k--')
    Cdist_grpd = paren(reshapeBOX(Cdist',subset),':',1);
    Ident_grpd = paren(reshapeBOX(Ident',subset),':',1);
    %Average, STD, n-1, t-value for segments larger than 100 bp
    ident.p(1) = mean(Ident_grpd(Cdist_grpd>100));
    ident.p(3) = numel(Ident_grpd(Cdist_grpd>100))-1;
    ident.p(2) = std(Ident_grpd(Cdist_grpd>100));
    ident.p(4) = (ident.p(1) - identAve)/(ident.p(2)/sqrt(ident.p(3)+1));
    %Average, STD, n-1, t-value for segments equal to or smaller than 100 bp
    ident.m(1) = mean(Ident_grpd(Cdist_grpd<=100));
    ident.m(3) = numel(Ident_grpd(Cdist_grpd<=100))-1;
    ident.m(2) = std(Ident_grpd(Cdist_grpd<=100));
    ident.m(4) = (ident.m(1) - identAve)/(ident.m(2)/sqrt(ident.m(3)+1));
    
    varargout{1} = ident;
end


j=0;
for i = subset
    j=j+1;
    subplot(4,4,[2:4 6:8 10:12])
    scatter(Cdist{i},Ident{i},24,clr(j,:),'filled','MarkerFaceAlpha',.7); 
    hold on;
end
 
ylim([ylimit]); xlim([xlimit]);
grid on; box on
set(gca,'XTickLabel',[],'YTickLabel',[]);
A = exist('identAve');
if A~=0
    if identAve>1
        identAve = identAve/100;
    end
    hold on; plot([xlimit],[identAve identAve],'k--')
end

[cnts, bins, cnts_idx] = histcounts(paren(reshapeBOX(Cdist',subset),':',1),'BinLimits',xlimit);
errorbar((bins(1:end-1)+bins(2)/2)',accumarray(cnts_idx, paren(reshapeBOX(Ident',subset),':',1),[], @mean),...
    accumarray(cnts_idx, paren(reshapeBOX(Ident',subset),':',1),[], @std),'ro','linewidth',1.2);

