%plot2.m % Load all CNPSummaries before!!!
%clear all
close all

% Plot all clusters of one example sample
savePath= "/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/Paper_kleinesPaper/plots_Isabel/";
printMe = "on";
basePath = "/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/Paper_kleinesPaper/allLists/";
recipSize = 4215607;
cmp = repmat([189 38 38; 0 76 255; 81 189 81; 232 163 23; 0 0 0;  235 235 84; 130 218 250; 247 134 247]/255, 70,1);
smpl = "BAns0620";
printName = "BATRO 6";
choseColor = cmp(3,:);

donor = "Batro";

if donor == "W23"
   accGenome = basePath + "accW232Bs166NCe.txt";
elseif donor == "Batro"
   accGenome = basePath + "accBatro2Bs166NCe.txt";
elseif donor == "Bval"
    accGenome = basePath + "accBval2Bs166NCe.txt";
end

fid = fopen(accGenome);
imp = textscan(fid, '%f %f %f', 'delimiter', '\t');
accStart = imp{1}; accEnd = imp{2};
fclose(fid); clear imp

if smpl~=""
    msk = contains([allCNPSummaries.Sample], smpl) & ~contains([allCNPSummaries.Sample], "No") ;
    adistStart = allCNPSummaries(msk).Adist(:,2);
    adistEnd = allCNPSummaries(msk).Adist(:,3);
end


figure(2); 
set(gcf, "Renderer", "painters", "Position", [750 250 400 400])
Lines = 6;

subPlotStepSize = round(recipSize / Lines,-1);

subplot(Lines, 1, 1); set(gca, "FontSize", 10); title(printName, "Interpreter", "none"); hold on;

for i = 1 : Lines
    
    
    subplot(Lines, 1, i); hold on;
    set(gca, "FontSize", 6)
    tmpStart = 0+(i-1)*subPlotStepSize; tmpEnd = subPlotStepSize+(i-1)*subPlotStepSize;
    plot([0+(i-1)*subPlotStepSize subPlotStepSize+(i-1)*subPlotStepSize], [0.5 0.5], "LineWidth", 15, "Color", [0.7 0.7 0.7 0.3])
    
    
    if smpl~=""
        %%% INTEGRATED SEGMENTS
        withIn = find(adistStart >= tmpStart & adistEnd <= tmpEnd);
            if ~isempty(withIn)
                p1 = plot([adistStart(withIn)'; adistEnd(withIn)'], [0.5 0.5], "LineWidth", 15, "Color", choseColor);
            end
       withOut = find(adistStart >= tmpStart & adistStart <= tmpEnd & adistEnd > tmpEnd);

       for j = 1 : length(withOut)
           adistStart = [adistStart; tmpEnd];
           adistEnd = [adistEnd; adistEnd(withOut(j))];
           plot([adistStart(withOut(j)) tmpEnd], [0.5 0.5], "LineWidth", 15, "Color",  choseColor)
       end
    end
    
   %%% ACCESSORY REGIONS
   withInAcc = find(accStart >= tmpStart & accEnd <= tmpEnd);
        if ~isempty(withInAcc)
            pAcc = plot([accStart(withInAcc) accEnd(withInAcc)], [0.5 0.5], "LineWidth", 15,  "Color", [0.5 0.5 0.5 0.8]);
        end
   withOutAcc = find(accStart >= tmpStart & accStart <= tmpEnd & accEnd > tmpEnd);
   
   for j = 1 : length(withOutAcc)
       accStart = [accStart; tmpEnd];
       accEnd = [accEnd; accEnd(withOutAcc(j))];
       plot([accStart(withOutAcc(j)) tmpEnd], [0.5 0.5], "LineWidth", 15, "Color", [0.5 0.5 0.5 0.8])
   end
   
   ylim([0.5-0.5 0.5+0.5])
   xlim([0+(i-1)*subPlotStepSize subPlotStepSize+(i-1)*subPlotStepSize])
   set(gca, "YTickLabel", [], "YTick", [])
   set(gca, "XTick", [tmpStart tmpEnd], "XTickLabel",[string(tmpStart); string(tmpEnd)] )
  
end 

% subplot(Lines, 1, Lines); 
% gca
legend([p1(1) pAcc(1)], "Det. Seg", "Not Core", "Location", "eastoutside")


xlabel("Position on recipient genome")

if printMe == "on"
    print(gcf,  '-painters', '-dpng', savePath + datestr(now, "yyyymmdd") + "_" + printName);
end