%plot2.m % Load all CNPSummaries before!!!

close all

% Plot all clusters of one example sample
basePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/1_kleinesPaper/allLists/";
printMe = "off";

recipSize = 4215607;


donors     = ["B. spiz.","B. val.","B. atro.","G. thermog."];
accGenome  = basePath + "acc/" + ["accW23_2NCe.txt","accBval_2NCe.txt","accBatro_2NCe.txt","accGeo_2NCe.txt"];
collectAcc = struct('donors',[],'accStart',[],'accEnd',[]);
ml         = basePath + "ml/" + ["mlW232Bs166NCe_v1.mat","mlBval2Bs166NCe_v1.mat","mlBatro2Bs166NCe_v1.mat","mlGeo2Bs166NCe_v1.mat"];
collectml  = struct('donors',[],'masterlist',[]);

hists = []; histLegends = []; histmax = []; histBinWidth = 0.002;
histBinEdges = -40*histBinWidth-0.5*histBinWidth:histBinWidth:40*histBinWidth+0.5*histBinWidth;
%cmp = repmat([189 38 38; 0 76 255; 81 189 81; 232 163 23; 0 0 0;  235 235 84; 130 218 250; 247 134 247]/255, 70,1);

DonorColor = [51 34 136; 136 34 85; 68 170 153; 205 102 119]/265;
gray = [221 221 221]/265;

%% load data

for i=1:numel(donors)
    collectAcc(i).donors = donors(i);
    
    fid = fopen(accGenome(i));
    imp = textscan(fid, '%f %f %f', 'delimiter', '\t');
    accStart = imp{1}; accEnd = imp{2};
    
    
    collectAcc(i).accStart = accStart;
    collectAcc(i).accEnd   = accEnd;
    
    % more interesting stuff
    collectAcc(i).accL   = accEnd - accStart + 1;
    collectAcc(i).coreL  = [accStart(1); accStart(2:end) - accEnd(1:end-1); recipSize - accEnd(end)];
    
    collectAcc(i).fracAcc = sum(collectAcc(i).accL)/recipSize;
    fclose(fid); clear imp; clear accStart accEnd
    
    % % now look at masterlist
    
    collectml(i).donors = donors(i);
    load(ml(i));
    collectml(i).masterlist = mlSummary.masterlist;
    clear mlSummary
end




savePath = "/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/Paper_kleinesPaper/plots_Isabel/";

%% start with plotting
for s=1:numel(donors)
    
    figure(s);
    set(gcf, "Renderer", "painters", "Position", [750 250 400 400])
    subPlots = 6;
    
    subPlotStepSize = round(recipSize / subPlots,-1);
    subplot(subPlots, 1, 1); set(gca, "FontSize", 10); hold on;
    
    accStart = collectAcc(s).accStart;
    accEnd   = collectAcc(s).accEnd;

    
    for i = 1 : subPlots
        
        if i==1
            title("accessory genome of " + donors(s))
        end
        subplot(subPlots, 1, i); hold on;
        set(gca, "FontSize", 6)
        tmpStart = 0+(i-1)*subPlotStepSize; tmpEnd = subPlotStepSize+(i-1)*subPlotStepSize;
        plot([0+(i-1)*subPlotStepSize subPlotStepSize+(i-1)*subPlotStepSize], [0.5 0.5], "LineWidth", 15, "Color", [0.9 0.9 0.9 0.3])
        
        
        %%% ACCESSORY REGIONS
        withInAcc = find(accStart >= tmpStart & accEnd <= tmpEnd);
        if ~isempty(withInAcc)
            pAcc = plot([accStart(withInAcc) accEnd(withInAcc)], [0.5 0.5], "LineWidth", 15,  "Color", DonorColor(s,:));
        end
        withOutAcc = find(accStart >= tmpStart & accStart <= tmpEnd & accEnd > tmpEnd);
        
        for j = 1 : length(withOutAcc)
            accStart = [accStart; tmpEnd];
            accEnd = [accEnd; accEnd(withOutAcc(j))];
            plot([accStart(withOutAcc(j)) tmpEnd], [0.5 0.5], "LineWidth", 15, "Color", DonorColor(s,:));
        end
        
        ylim([0.5-0.5 0.5+0.5])
        xlim([0+(i-1)*subPlotStepSize subPlotStepSize+(i-1)*subPlotStepSize])
        set(gca, "YTickLabel", [], "YTick", [])
        set(gca, "XTick", [tmpStart tmpEnd], "XTickLabel",[string(tmpStart); string(tmpEnd)] )
        
        %     set(gca, "Visible", "off")
    end
    
    xlabel("Position on recipient genome")
    if printMe == "on"
        print(figure(s),  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_" + donors(s) + "_Acc");
    end
    
end

%% analysing more stuff

% distributions of core and acc piece lengths
figure(10); hold on; box on;
set(gcf, 'Position', [750 250 500 500], 'Renderer', 'painters')
histBinWidth = 1000;
histBinEdges = 0:histBinWidth:35*histBinWidth;

for i=1:numel(donors)
subplot(4,1,i); hold on;
set(gca, 'FontSize', 9,'YScale','log');
h = histogram(collectAcc(i).accL, 'BinEdges', histBinEdges,...
    'Normalization','probability','FaceColor','none','EdgeColor',DonorColor(i,:),'LineWidth',1.5,...
    'DisplayName',donors(i));
l = plot([median([collectAcc(i).accL]) median([collectAcc(i).accL])],[0.001 1.1],'k--','LineWidth',2,...
    'DisplayName','median');
legend
ylim([0 1.1])

if i == numel(donors)
xlabel('length of accessory parts')
end


end

figure(11); hold on; box on;
set(gcf, 'Position', [0 0 500 500], 'Renderer', 'painters')

histBinWidth = 1000;
histBinEdges = 0:histBinWidth:45*histBinWidth;

for i=1:numel(donors)
    subplot(4,1,i); hold on;
    set(gca, 'FontSize', 9,'YScale','log');
    h = histogram(collectAcc(i).coreL, 'BinEdges', histBinEdges,...
        'Normalization','probability','FaceColor','none','EdgeColor',DonorColor(i,:),'LineWidth',1.5,...
        'DisplayName',donors(i));
    l = plot([median([collectAcc(i).coreL]) median([collectAcc(i).coreL])],[0.001 1.1],'k--','LineWidth',2,...
        'DisplayName','median');
    legend
    ylim([0 1.1])
    xlim([0 30000])
    
    if i == numel(donors)
        xlabel('length of core parts')
    end
end

%% analyse stuff with the masterlist

identity = [];
for i=1:numel(donors)
   identity(i).donors = collectml(i).donors;
   
   identity(i).ident  =  1 - length(collectml(i).masterlist)/ (recipSize - sum([collectAcc(i).accL]));
    
end


%% printing the plots
if printMe == "on"
    print(figure(10), '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_" + donor + "_LengthDistAcc");
    print(figure(11), '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_" + donor + "_LengthDistCore");
end