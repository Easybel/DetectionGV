close all

Wc10Mask = contains([allCNPSummaries.Sample], "Wns") & contains([allCNPSummaries.Sample], "10")& ~contains([allCNPSummaries.Sample], "No");
Vc10Mask = contains([allCNPSummaries.Sample], "Vns") & contains([allCNPSummaries.Sample], "10") & ~contains([allCNPSummaries.Sample], "No");
Bc10Mask = contains([allCNPSummaries.Sample], "BAns") & contains([allCNPSummaries.Sample], "10");
Noc10Mask = contains([allCNPSummaries.Sample], "No") & contains([allCNPSummaries.Sample], "10");
recipSize = 4215607;
cmp = repmat([189 38 38; 0 76 255; 81 189 81; 232 163 23; 0 0 0;  235 235 84; 130 218 250; 247 134 247]/255, 70,1);


%whichAcc = "BlastN";
whichAcc = "ours";

% Average number of segments per replicate:
mean(cellfun(@(x) length(x), {allCNPSummaries(Wc10Mask).Cdist}))
std(cellfun(@(x) length(x), {allCNPSummaries(Wc10Mask).Cdist}))

mean(cellfun(@(x) length(x), {allCNPSummaries(Vc10Mask).Cdist}))
std(cellfun(@(x) length(x), {allCNPSummaries(Vc10Mask).Cdist}))

mean(cellfun(@(x) length(x), {allCNPSummaries(Bc10Mask).Cdist}))
std(cellfun(@(x) length(x), {allCNPSummaries(Bc10Mask).Cdist}))

% Average length of segments
adistW = [];
for i = 1 : length(nonzeros(Wc10Mask))
    tmpCNP = allCNPSummaries(Wc10Mask);
    adistW = [adistW; tmpCNP(i).Adist(:,1)];
end
mean(adistW)

adistV = [];
for i = 1 : length(nonzeros(Vc10Mask))
    tmpCNP = allCNPSummaries(Vc10Mask);
    adistV = [adistV; tmpCNP(i).Adist(:,1)];
end
mean(adistV)

adistA = [];
for i = 1 : length(nonzeros(Bc10Mask))
    tmpCNP = allCNPSummaries(Bc10Mask);
    if ~isempty(tmpCNP(i).Adist)
        adistA = [adistA; tmpCNP(i).Adist(:,1)];
    end
end
mean(adistA)


%%  Percentage of genome replace
savePath  = "/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/Paper_kleinesPaper/plots_Isabel/";

Wc20Mask = contains([allCNPSummaries.Sample], "Wns") & contains([allCNPSummaries.Sample], "20")& ~contains([allCNPSummaries.Sample], "No");
Vc20Mask = contains([allCNPSummaries.Sample], "Vns") & contains([allCNPSummaries.Sample], "20") & ~contains([allCNPSummaries.Sample], "No");
Bc20Mask = contains([allCNPSummaries.Sample], "BAns") & contains([allCNPSummaries.Sample], "20");
Noc20Mask = contains([allCNPSummaries.Sample], "No") & contains([allCNPSummaries.Sample], "20");

transW10 = [allCNPSummaries(Wc10Mask).Transfer]/0.91;
transW20 = [allCNPSummaries(Wc20Mask).Transfer]/0.91;
transV10 = [allCNPSummaries(Vc10Mask).Transfer]/0.84;
transV20 = [allCNPSummaries(Vc20Mask).Transfer]/0.84;
transB10 = [allCNPSummaries(Bc10Mask).Transfer]/0.69;
transB20 = [allCNPSummaries(Bc20Mask).Transfer]/0.69;

% Our core/acc
if whichAcc == "ours"
    transW10 = [allCNPSummaries(Wc10Mask).Transfer]/0.87;
    transW20 = [allCNPSummaries(Wc20Mask).Transfer]/0.87;
    transV10 = [allCNPSummaries(Vc10Mask).Transfer]/0.86;
    transV20 = [allCNPSummaries(Vc20Mask).Transfer]/0.86;
    transB10 = [allCNPSummaries(Bc10Mask).Transfer]/0.75;
    transB20 = [allCNPSummaries(Bc20Mask).Transfer]/0.75;
end

ursprung = zeros(1, 1);

figure(3); hold on; box on; title("\itB. spizizenii")
set(gcf, "Renderer", "painters", "Position", [25 570 476 242]);
set(gca, "FontSize", 12)

% allCombW = combvec(ursprung, transW10, transW20);
m = 0;
ft1 = fittype({'x'}); % no intercept f(a,x) = ax
for j = 1 : length(transW10)
    for k = 1: length(transW20)
        m = m+1;
        dataPoints = [0 transW10(j) transW20(k)];
        Fit1 = fit([0; 20; 40], dataPoints', ft1);
        
        F1(m) = Fit1.a;
    end
end

xPlot = 0: 1 : 60;
yF1 = mean(F1(:))*xPlot;
yF1plus = (mean(F1(:))+std(F1(:)))*xPlot;
yF1minus = (mean(F1(:))-std(F1(:)))*xPlot;

pFit1 = plot(xPlot, yF1, "-", "Color", [0.7 0.7 0.9], "LineWidth", 1.4);
plot(xPlot, yF1plus, "--", "Color", [0.7 0.7 0.9], "LineWidth", 1.4);
plot(xPlot, yF1minus, "--", "Color", [0.7 0.7 0.9], "LineWidth", 1.4);

if whichAcc == "BlastN"
    plot(repmat(20, 1, length(nonzeros(Wc10Mask))), [allCNPSummaries(Wc10Mask).Transfer]/0.91, "x", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8]); %, "Color", [0.7 0.7 0.7 0.8])
    plot(20, mean([allCNPSummaries(Wc10Mask).Transfer]/0.91), "x", "LineWidth", 1.4, "Color", [0 0 0])
    plot(repmat(40, 1, length(nonzeros(Wc20Mask))), [allCNPSummaries(Wc20Mask).Transfer]/0.91, "x", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8])
    plot(40, mean([allCNPSummaries(Wc20Mask).Transfer]/0.91), "x", "LineWidth", 1.4, "Color", [0 0 0])
    
elseif whichAcc == "ours"
    plot(repmat(20, 1, length(nonzeros(Wc10Mask))), [allCNPSummaries(Wc10Mask).Transfer]/0.91, "x", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8]); %, "Color", [0.7 0.7 0.7 0.8])
    plot(20, mean([allCNPSummaries(Wc10Mask).Transfer]/0.87), "x", "LineWidth", 1.4, "Color", [0 0 0])
    plot(repmat(40, 1, length(nonzeros(Wc20Mask))), [allCNPSummaries(Wc20Mask).Transfer]/0.91, "x", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8])
    plot(40, mean([allCNPSummaries(Wc20Mask).Transfer]/0.87), "x", "LineWidth", 1.4, "Color", [0 0 0])
end

xlim([0 60]); xlabel("DNA uptake [h]"); ylabel("Genome replaced [%]")
legend(pFit1, "Rate: " + num2str(mean(F1(:)), "%.2f") + " \pm " +num2str(std(F1(:)), "%.2f")+ " % per hour ", "Location", "northwest")
ylim([0 15])
print(gcf,  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_W_genRepl");




figure(4); hold on; box on; title("\itB. vallismortis")
set(gcf, "Renderer", "painters", "Position", [25 570 476 242]);
set(gca, "FontSize", 12)

%allCombV = combvec(ursprung, transV10, transV20);
m = 0;
ft1 = fittype({'x'}); % no intercept f(a,x) = ax
for j = 1 : length(transV10)
    for k = 1:length(transV20)
        m = m + 1;
        dataPoints = [0 transV10(j) transV20(k)];
        Fit1 = fit([0; 20; 40], dataPoints', ft1);
        F2(m) = Fit1.a;
    end
end

xPlot = 0: 1 : 60;
yF2 = mean(F2(:))*xPlot;
yF2plus = (mean(F2(:))+std(F2(:)))*xPlot;
yF2minus = (mean(F2(:))-std(F2(:)))*xPlot;

pFit2 = plot(xPlot, yF2, "-", "Color", [0.7 0.7 0.9], "LineWidth", 1.4);
plot(xPlot, yF2plus, "--", "Color", [0.7 0.7 0.9], "LineWidth", 1.4);
plot(xPlot, yF2minus, "--", "Color", [0.7 0.7 0.9], "LineWidth", 1.4);

if whichAcc == "BlastN"
    plot(repmat(20, 1, length(nonzeros(Vc10Mask))), [allCNPSummaries(Vc10Mask).Transfer]/0.84, "kx", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8])
    plot(20, mean([allCNPSummaries(Vc10Mask).Transfer]/0.84), "x", "LineWidth", 1.4, "Color", [0 0 0])
    plot(repmat(40, 1, length(nonzeros(Vc20Mask))), [allCNPSummaries(Vc20Mask).Transfer]/0.84, "kx", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8])
    plot(40, mean([allCNPSummaries(Vc20Mask).Transfer]/0.84), "x", "LineWidth", 1.4, "Color", [0 0 0])
    
elseif whichAcc == "ours"
    plot(repmat(20, 1, length(nonzeros(Vc10Mask))), [allCNPSummaries(Vc10Mask).Transfer]/0.84, "kx", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8])
    plot(20, mean([allCNPSummaries(Vc10Mask).Transfer]/0.86), "x", "LineWidth", 1.4, "Color", [0 0 0])
    plot(repmat(40, 1, length(nonzeros(Vc20Mask))), [allCNPSummaries(Vc20Mask).Transfer]/0.84, "kx", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8])
    plot(40, mean([allCNPSummaries(Vc20Mask).Transfer]/0.86), "x", "LineWidth", 1.4, "Color", [0 0 0])
end

xlim([0 60]); xlabel("DNA uptake [h]"); ylabel("Genome replaced [%]")

legend(pFit2, "Rate: " + num2str(mean(F2(:)), "%.2f") + " \pm " +num2str(std(F2(:)), "%.2f")+ " % per hour ", "Location", "northwest")
ylim([0 9])
print(gcf,  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_V_genRepl");



figure(5); hold on; box on; title("\itB. atrophaeus")
set(gcf, "Renderer", "painters", "Position", [25 570 476 242]);
set(gca, "FontSize", 12)

m = 0;
ft1 = fittype({'x'}); % no intercept f(a,x) = ax
for j = 1 : length(transB10)
    for k = 1: length(transB20)
        m = m + 1;
        dataPoints = [0 transB10(j) transB20(k)];
        Fit1 = fit([0; 20; 40], dataPoints', ft1);
        F3(m) = Fit1.a;
    end
end

xPlot = 0: 1 : 60;
yF3 = mean(F3(:))*xPlot;
yF3plus = (mean(F3(:))+std(F3(:)))*xPlot;
yF3minus = (mean(F3(:))-std(F3(:)))*xPlot;

pFit3 = plot(xPlot, yF3, "-", "Color", [0.7 0.7 0.9], "LineWidth", 1.4);
plot(xPlot, yF3plus, "--", "Color", [0.7 0.7 0.9], "LineWidth", 1.4);
plot(xPlot, yF3minus, "--", "Color", [0.7 0.7 0.9], "LineWidth", 1.4);

if whichAcc == "BlastN"
    plot(repmat(20, 1, length(nonzeros(Bc10Mask))), [allCNPSummaries(Bc10Mask).Transfer]/0.69, "x", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8])
    plot(20, mean([allCNPSummaries(Bc10Mask).Transfer]/0.69), "x", "LineWidth", 1.4, "Color", [0 0 0])
    plot(repmat(40, 1, length(nonzeros(Bc20Mask))), [allCNPSummaries(Bc20Mask).Transfer]/0.69, "kx", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8])
    plot(40, mean([allCNPSummaries(Bc20Mask).Transfer]/0.69), "x", "LineWidth", 1.4, "Color", [0 0 0])
    
elseif whichAcc == "ours"
    plot(repmat(20, 1, length(nonzeros(Bc10Mask))), [allCNPSummaries(Bc10Mask).Transfer]/0.69, "x", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8])
    plot(20, mean([allCNPSummaries(Bc10Mask).Transfer]/0.75), "x", "LineWidth", 1.4, "Color", [0 0 0])
    plot(repmat(40, 1, length(nonzeros(Bc20Mask))), [allCNPSummaries(Bc20Mask).Transfer]/0.69, "kx", "LineWidth", 1.4, "Color", [0.7 0.7 0.7 0.8])
    plot(40, mean([allCNPSummaries(Bc20Mask).Transfer]/0.75), "x", "LineWidth", 1.4, "Color", [0 0 0])
end

xlim([0 60]); xlabel("DNA uptake [h]"); ylabel("Genome replaced [%]")

legend(pFit3, "Rate: " + num2str(mean(F3(:)), "%.3f") + " \pm " +num2str(std(F3(:)), "%.3f")+ " % per hour ", "Location", "northwest")
ylim([0 0.5])
print(gcf,  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_B_genRepl");



%%
bpAcc      = [538402 604898 1071911 3738778];
bpCore     = recipSize - [538402 604898 1071911 3738778];
% SNPs that are differnt in core genome = mlSNPs and that are the same
% identSNPs
mlSNPs     = [249596 300952 424514 94128];
identSNPs  = bpCore - mlSNPs;

div        = 100 * mlSNPs ./ bpCore;
spec_ident = 100 - div;

rate      = [mean(F1) mean(F2) mean(F3) 0];
rateerror = [std(F1) std(F2) std(F3) 0];
rate_bp   = rate.*bpCore/100;
rateerror_bp = rateerror .* bpCore /100 ;

clr = [0.2 0.8 0.1; 0.4 0.2 0.9; 0.9 0.5 0.2; 0.1 0.1 0.1];

% % %
figure(6); hold on; box on; title("")
set(gcf, "Renderer", "painters", "Position", [25 570 500 450]);
set(gca, "FontSize", 12)

for i=1:3
p(i) = plot(spec_ident(i) * bpCore(i) /100, rate(i), "x", "LineWidth", 1.4, "Color",cmp(i,:));
errorbar(spec_ident(i) * bpCore(i)/100, rate(i),rateerror(i), "LineWidth", 1.4, "Color", cmp(i,:));
end
p(4)  = plot([spec_ident(4) * bpCore(4)/100 spec_ident(4) * bpCore(4)/100], [10e-18 100], '--k');

ft2 = fittype({'x', '1'});
Fit1 = fit(spec_ident(1:3)'.*bpCore(1:3)'/100, log(rate(1:3)'), ft2);
Fit11 = fitlm()
PRate(1) = Fit1.a;
PRate(2) = Fit1.b;

xPlot = 0 :100: 4000000;
yfit = exp(PRate(1) * xPlot + PRate(2));
plot(xPlot, yfit, 'k-.');
set(gca, "YScale", "log")

xlim([2500000 3600000])
ylim([5e-4 1])
legend([p(1) p(2) p(3) p(4)], ["\itB. spizizenii", "\itB. vallismortis", "\itB. atrophaeus","\itGeo. thermo."],'Location','NorthWest')

ylabel("Transfer rate in\newline core genome [% per hour]")
xlabel("Average Identity * core genome \newline [bp]")
print(figure(6),  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_rateVsIdentinbpZoom");
figure(6)
xlim([6 15])
ylim([5e-4 1])
print(figure(6),  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_rateVsDivZoom");

% fitgüte
fit_goodness_identCore = sum((exp(PRate(1) * (spec_ident(1:3)'.*bpCore(1:3)'/100) + PRate(2)) - rate(1:3)').^2);    

% % % %
figure(61); hold on; box on; title("")
set(gcf, "Renderer", "painters", "Position", [25 570 500 450]);
set(gca, "FontSize", 12)

for i=1:3
p(i) = plot(spec_ident(i), rate(i), "x", "LineWidth", 1.4, "Color",cmp(i,:));
errorbar(spec_ident(i), rate(i),rateerror(i), "LineWidth", 1.4, "Color", cmp(i,:));
end
p(4)  = plot([spec_ident(4) spec_ident(4)], [10e-18 100], '--k');

ft2 = fittype({'x', '1'});
Fit1 = fit(spec_ident(1:3)', log(rate(1:3)'), ft2);
PRate(1) = Fit1.a;
PRate(2) = Fit1.b;
xPlot = 0 :0.1: 100;
yfit = exp(PRate(1) * xPlot + PRate(2));
plot(xPlot, yfit, 'k-.');
set(gca, "YScale", "log")

xlim([85 95])
ylim([5e-4 1])
legend([p(1) p(2) p(3) p(4)], ["\itB. spizizenii", "\itB. vallismortis", "\itB. atrophaeus","\itGeo. thermo."],'Location','NorthWest')

ylabel("Transfer rate in\newline core genome [% per hour]")
xlabel("Average Identity of core genome \newline [bp]")
print(figure(61),  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_rateVsIdentinbpZoom");
figure(61)
xlim([6 15])
ylim([5e-4 1])
print(figure(61),  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_rateVsDivZoom");

% fitgüte

fit_goodness_ident = sum((exp(PRate(1) * (spec_ident(1:3)') + PRate(2)) - rate(1:3)').^2);   

%%
figure(7); hold on; box on; title("")
set(gcf, "Renderer", "painters", "Position", [25 570 500 450]);
set(gca, "FontSize", 12)

for i=1:3
p(i) = plot(mlSNPs(i), rate(i), "x", "LineWidth", 1.4, "Color",clr(i,:));
errorbar(mlSNPs(i), rate(i),rateerror(i), "Color", clr(i,:));
end
p(4)  = plot([mlSNPs(4) mlSNPs(4)], [10e-7 1e10], '--k');

ft2 = fittype({'x', '1'});
Fit1 = fit(mlSNPs(1:3)', log(rate(1:3)'), ft2);
PRate(1) = Fit1.a;
PRate(2) = Fit1.b;
xPlot = 0 :1000: 4000000;
yfit = exp(PRate(1) * xPlot + PRate(2));
plot(xPlot, yfit, 'k-.');
set(gca, "YScale", "log")

xlim([90000 450000])
ylim([1e-5 1.1])

ylabel("Transfer Rate\newline[% per hour]")
xlabel("Length of core genome")

legend([p(1) p(2) p(3) p(4)], ["\itB. spizizenii", "\itB. vallismortis", "\itB. atrophaeus","\itGeo. thermoglucosidasius"],'Location','NorthWest')

print(gcf,  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_CoreLengthVsDiv");

%% look at identites
idx_Wc20Mask = find(Wc20Mask);
IdentW20_temp = []; LW20 = [];
idx_Vc20Mask = find(Vc20Mask);
IdentV20_temp = []; LV20 = [];
idx_Bc20Mask = find(Bc20Mask);
IdentB20_temp = []; LB20 = [];

for i=1:numel(idx_Wc20Mask)
   IdentW20_temp =  [IdentW20_temp; allCNPSummaries(idx_Wc20Mask(i)).Ident];
   LW20     =  [LW20; allCNPSummaries(idx_Wc20Mask(i)).Adist(:,1)]
end

for i=1:numel(idx_Vc20Mask)
   IdentV20_temp =  [IdentV20_temp; allCNPSummaries(idx_Vc20Mask(i)).Ident];
   LV20     =  [LV20; allCNPSummaries(idx_Vc20Mask(i)).Adist(:,1)]
end

for i=1:numel(idx_Bc20Mask)
   IdentB20_temp =  [IdentB20_temp; allCNPSummaries(idx_Bc20Mask(i)).Ident];
   if ~isempty(allCNPSummaries(idx_Bc20Mask(i)).Adist)
   LB20     =  [LB20; allCNPSummaries(idx_Bc20Mask(i)).Adist(:,1)]
   end
end

LW20_mask = LW20 > 100;
LV20_mask = LV20 > 100;
LB20_mask = LB20 > 180;

IdentW20 = IdentW20_temp(LW20_mask);
IdentV20 = IdentV20_temp(LV20_mask);
IdentB20 = IdentB20_temp(LB20_mask);

figure(20); hold on; box on;
set(gcf, 'Position', [0 0 500 500], 'Renderer', 'painters')

histBinWidth = 0.01;
histBinEdges = 0:histBinWidth:100*histBinWidth;

subplot(3,1,1); hold on;
%set(gca, 'FontSize', 9,'YScale','log');
h1 = histogram(IdentW20, 'BinEdges', histBinEdges,...
    'Normalization','probability','FaceColor','none','EdgeColor',cmp(1,:),'LineWidth',1.5);
l(1) = plot([mean(IdentW20) mean(IdentW20)],[0 1.1],'LineStyle','--','Color',cmp(1,:),'LineWidth',2);
ylim([0 0.4])
xlim([0.75 1.01])
legend([h1 l(1)], "BSPIZ cy20",...
    "Mean = " + num2str(mean(IdentW20), "%.4f") + "\newlineStd = " + num2str(std(IdentW20), "%.4f"),...
    'Location','NorthWest')

subplot(3,1,2); hold on;
h2 = histogram(IdentV20, 'BinEdges', histBinEdges,...
    'Normalization','probability','FaceColor','none','EdgeColor',cmp(2,:),'LineWidth',1.5);
l(2) = plot([mean(IdentV20) mean(IdentV20)],[0 1.1],'LineStyle','--','Color',cmp(2,:),'LineWidth',2);
ylim([0 0.4])
xlim([0.75 1.01])
legend([h2 l(2)],"BVAL cy20",...
    "Mean = " + num2str(mean(IdentV20), "%.4f") + "\newlineStd = " + num2str(std(IdentV20), "%.4f"),...
    'Location','NorthWest')

subplot(3,1,3); hold on;
h3 = histogram(IdentB20, 'BinEdges', histBinEdges,...
    'Normalization','probability','FaceColor','none','EdgeColor',cmp(3,:),'LineWidth',1.5);
l(3) = plot([mean(IdentB20) mean(IdentB20)],[0 1.1],'LineStyle','--','Color',cmp(3,:),'LineWidth',2);
legend([h3 l(3)],"BATRO cy20",...
    "Mean = " + num2str(mean(IdentB20), "%.4f") + "\newlineStd = " + num2str(std(IdentB20), "%.4f"),...
    'Location','NorthWest')
ylim([0 0.4])
xlim([0.75 1.01])


ylabel('frequency')
xlabel('Identity of segements > 100 bp')
print(figure(20),  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_Ident_compare");

%% do a bootstrap - yippieh!!
N = 10000;
[W23_bootstat,W23_sample] = bootstrp(N,@(x)[mean(x) std(x)],IdentW20);
[Bval_bootstat,Bval_sample] = bootstrp(N,@(x)[mean(x) std(x)],IdentV20);
[Batro_bootstat,Batro_sample] = bootstrp(N,@(x)[mean(x) std(x)],IdentB20);

%% plot the results
%   -- the botstrap mean and std distribution

histBinW(1) = 0.001;
histBinW(2) = 0.001;

for i=1:2
    
    figure(40 + i); set(gcf, 'Position', [550 350 700 420], 'Renderer', 'painters');
    hold on; box on;
    set(gca, 'FontSize', 9);
    

    histBinEdges = 0:histBinW(i):1000*histBinW(i);
    
    h_bootW23 = histogram(W23_bootstat(:,i), 'BinEdges', histBinEdges,...
        'Normalization', 'probability','FaceColor','none','EdgeColor',cmp(1,:),'LineWidth',1.5);
    h_bootBval = histogram(Bval_bootstat(:,i), 'BinEdges', histBinEdges,...
        'Normalization', 'probability','FaceColor','none','EdgeColor',cmp(2,:),'LineWidth',1.5);
%     h_bootBatro = histogram(Batro_bootstat(:,i), 'BinEdges', histBinEdges,...
%         'Normalization', 'probability','FaceColor',cmp(3,:),'FaceAlpha', 0.4, 'EdgeColor',cmp(3,:),'LineWidth',1.5);
    
    MbootCW23  = mean(W23_bootstat(:,i));
    MbootCBval = mean(Bval_bootstat(:,i));
    MbootCBatro  = mean(Batro_bootstat(:,i));
    
    SbootCW23  =  std(W23_bootstat(:,i));
    SbootCBval =  std(Bval_bootstat(:,i));
%     SbootCBatro  =  std(Batro_bootstat(:,i));
    
    p_MbootW23 = plot([MbootCW23 MbootCW23], [0 2], "--",...
        "LineWidth", 1.2, "Color", cmp(1,:));
    p_MbootBval = plot([MbootCBval MbootCBval], [0 2], "--",...
        "LineWidth", 1.2, "Color", cmp(2,:));
%     p_MbootBatro = plot([MbootCBatro MbootCBatro], [0 2], "--",...
%         "LineWidth", 1.2, "Color", [0.2 0.2 0.2 0.5]);
    
    legend([h_bootW23 h_bootBval p_MbootW23 p_MbootBval],"BSPIZ cy 20", "BVAL cy20",...
        "Mean = " + num2str(MbootCW23, "%.4f") + "\newlineStd = " + num2str(SbootCW23, "%.4f"),...
        "Mean = " + num2str(MbootCBval, "%.4f") + "\newlineStd = " + num2str(SbootCBval, "%.4f"),...
        'Location','NorthWest')
    
    ylabel("frequency");
    
    if i==1
        xlim([0.9 0.95])
        ylim([0 0.6])
        xlabel("mean of the bootstrap samples");
        title("Bootstrapping with Identity: Sampling distribution of the sample mean")
        
    elseif i==2
        
        xlim([0 0.04])
        ylim([0 0.6])
        xlabel("standard deviation of the bootstrap samples");
        title("Bootstrapping with Identity: Sampling distribution of the sample std")
    end
    
end

%% find the confidence intervals of the bootstrap
alpha = 0.05;
z = 1.96;
interval = (N * alpha)/2;

W23_bootRes.meanMean   = mean(W23_bootstat(:,1));
W23_bootRes.stdMean    =  std(W23_bootstat(:,1));
W23_bootRes.meanStd    = mean(W23_bootstat(:,2));
W23_bootRes.stdStd     =  std(W23_bootstat(:,2));
sortW23Mean = sort(W23_bootstat(:,1));
W23_bootRes.MeanCIRank  = sort([sortW23Mean(interval) sortW23Mean(N-interval)]);
W23_bootRes.MeanCIwMean = [W23_bootRes.meanMean - z*W23_bootRes.stdMean W23_bootRes.meanMean + z*W23_bootRes.stdMean];
sortW23Std = sort(W23_bootstat(:,2));
W23_bootRes.StdCIRank  = sort([sortW23Std(interval) sortW23Std(N-interval)]);
W23_bootRes.StdCIwMean = [W23_bootRes.meanStd - z*W23_bootRes.stdStd W23_bootRes.meanStd + z*W23_bootRes.stdStd];

Bval_bootRes.meanMean   = mean(Bval_bootstat(:,1));
Bval_bootRes.stdMean    =  std(Bval_bootstat(:,1));
Bval_bootRes.meanStd    = mean(Bval_bootstat(:,2));
Bval_bootRes.stdStd     =  std(Bval_bootstat(:,2));
sortBvalMean = sort(Bval_bootstat(:,1));
Bval_bootRes.MeanCIRank  = sort([sortBvalMean(interval) sortBvalMean(N-interval)]);
Bval_bootRes.MeanCIwMean = [Bval_bootRes.meanMean - z*Bval_bootRes.stdMean Bval_bootRes.meanMean + z*Bval_bootRes.stdMean];
sortBvalStd = sort(Bval_bootstat(:,2));
Bval_bootRes.StdCIRank  = sort([sortBvalStd(interval) sortBvalStd(N-interval)]);
Bval_bootRes.StdCIwMean = [Bval_bootRes.meanStd - z*Bval_bootRes.stdStd Bval_bootRes.meanStd + z*Bval_bootRes.stdStd];

% Batro_bootRes.meanMean   = mean(Batro_bootstat(:,1));
% Batro_bootRes.stdMean    =  std(Batro_bootstat(:,1));
% Batro_bootRes.meanStd    = mean(Batro_bootstat(:,2));
% Batro_bootRes.stdStd     =  std(Batro_bootstat(:,2));
% sortBatroMean = sort(Batro_bootstat(:,1));
% Batro_bootRes.MeanCIRank  = sort([sortBatroMean(interval) sortBatroMean(N-interval)]);
% Batro_bootRes.MeanCIwMean = [Batro_bootRes.meanMean - z*Batro_bootRes.stdMean Batro_bootRes.meanMean + z*Batro_bootRes.stdMean];
% sortBatroStd = sort(Batro_bootstat(:,2));
% Batro_bootRes.StdCIRank  = sort([sortBatroStd(interval) sortBatroStd(N-interval)]);
% Batro_bootRes.StdCIwMean = [Batro_bootRes.meanStd - z*Batro_bootRes.stdStd Batro_bootRes.meanStd + z*Batro_bootRes.stdStd];

%%
% order: W23, Bval, Batro
height = [0.52 0.54 0.56];


figure(41);
hold on;

p_MbootW23 = plot(W23_bootRes.MeanCIwMean, [height(1) height(1)], '-',...
    "LineWidth", 1.2, "Color", cmp(1,:),'DisplayName','95%-CI');
plot(W23_bootRes.MeanCIwMean, [height(1) height(1)], "o",'Color', cmp(1,:),...
    'LineStyle','none','HandleVisibility','off');

p_MbootBval = plot(Bval_bootRes.MeanCIwMean, [height(2) height(2)], "-",...
    "LineWidth", 1.2, "Color", cmp(2,:),'DisplayName','95%-CI');
plot(Bval_bootRes.MeanCIwMean, [height(2) height(2)], "o",'Color', cmp(2,:),...
    'LineStyle','none','HandleVisibility','off');

% p_MbootBatro = plot(Batro_bootRes.MeanCIwMean, [height(3) height(3)], "-",...
%     "LineWidth", 1.2, "Color", cmp(3,:),'DisplayName','95%-CI');
% plot(Batro_bootRes.MeanCIwMean, [height(3) height(3)], "o",'Color',cmp(3,:),...
%     'LineStyle','none','HandleVisibility','off');


figure(42);
hold on;

p_MbootW23 = plot(W23_bootRes.StdCIwMean, [height(1) height(1)], '-',...
    "LineWidth", 1.2, "Color", cmp(1,:),'DisplayName','95%-CI');
plot(W23_bootRes.StdCIwMean, [height(1) height(1)], "o",'Color', cmp(1,:),...
    'LineStyle','none','HandleVisibility','off');

p_MbootBval = plot(Bval_bootRes.StdCIwMean, [height(2) height(2)], "-",...
    "LineWidth", 1.2, "Color", cmp(2,:),'DisplayName','95%-CI');
plot(Bval_bootRes.StdCIwMean, [height(2) height(2)], "o",'Color', cmp(2,:),...
    'LineStyle','none','HandleVisibility','off');
% 
% p_MbootBatro = plot(Batro_bootRes.StdCIwMean, [height(3) height(3)], "-",...
%     "LineWidth", 1.2, "Color", cmp(3,:),'DisplayName','95%-CI');
% plot(Batro_bootRes.StdCIwMean, [height(3) height(3)], "o",'Color', cmp(3,:),...
%     'LineStyle','none','HandleVisibility','off');


%% print
print(figure(41),  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_bootstrap_Ident_mean");

print(figure(42),  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_bootstrap_Ident_std");

