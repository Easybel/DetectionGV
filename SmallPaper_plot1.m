

Wc10Mask = contains([allSNPSummaries.Sample], "Wns") & contains([allSNPSummaries.Sample], "10")& ~contains([allSNPSummaries.Sample], "No");
Vc10Mask = contains([allSNPSummaries.Sample], "Vns") & contains([allSNPSummaries.Sample], "10") & ~contains([allSNPSummaries.Sample], "No"); 
Bc10Mask = contains([allSNPSummaries.Sample], "BAns") & contains([allSNPSummaries.Sample], "10");
Noc10Mask = contains([allCNPSummaries.Sample], "No") & contains([allCNPSummaries.Sample], "10");





clearvars boxPlot*
boxPlotCount = cellfun(@(x) length(x), {allSNPSummaries(Wc10Mask).IndvMutList}, "UniformOutput", true);
        noSNPmeanWc10 = round(mean(cellfun(@(x) length(x), {allSNPSummaries(Wc10Mask).IndvMutList}, "UniformOutput", true)));
boxPlotName = repmat("Spizizenii", 1, length(boxPlotCount));

boxPlotCount = [boxPlotCount cellfun(@(x) length(x), {allSNPSummaries(Vc10Mask).IndvMutList}, "UniformOutput", true)];
        noSNPmeanVc10 = round(mean(cellfun(@(x) length(x), {allSNPSummaries(Vc10Mask).IndvMutList}, "UniformOutput", true)));
boxPlotName = [boxPlotName repmat("Vallismortis", 1, length(boxPlotCount)-length(boxPlotName))];

boxPlotCount = [boxPlotCount cellfun(@(x) length(x), {allSNPSummaries(Bc10Mask).IndvMutList}, "UniformOutput", true)];
        noSNPmeanBc10 = mean(cellfun(@(x) length(x), {allSNPSummaries(Bc10Mask).IndvMutList}, "UniformOutput", true));
boxPlotName = [boxPlotName repmat("Atrophaeus", 1, length(boxPlotCount)-length(boxPlotName))]

boxPlotCount = [boxPlotCount cellfun(@(x) length(x), {allCNPSummaries(Noc10Mask).denovo}, "UniformOutput", true)];
        noSNPmeanNoc10 = mean(cellfun(@(x) length(x), {allCNPSummaries(Noc10Mask).denovo}, "UniformOutput", true));
boxPlotName = [boxPlotName repmat("Control", 1, length(boxPlotCount)-length(boxPlotName))]


figure(1); hold on;
set(gcf, "Renderer", "painters", "Position", [25 570 478 302]); 
set(gca, "FontSize", 10)
boxplot(boxPlotCount, boxPlotName)
ylabel("Number of detected SNPs")
title("After 10 cycles of the evolution experiment")
% set(gca, "YTickLabel", [0: 10000 : 40000])
ylim([1 1E5])
set(gca, "YScale", "log")
% set(gca, "XTickLabel"

pMean = plot([1 2 3 4], [noSNPmeanWc10 noSNPmeanVc10 noSNPmeanBc10 noSNPmeanNoc10], "rd", "MarkerSize", 3, "LineWidth", 3)

legend(pMean, "Mean")
print(gcf,  '-painters', '-dpng', SNPPath + datestr(now, "yyyymmdd") + "_numberSNPs_logScale");



