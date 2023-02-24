%% load the SNP/CNP Summaries

clearvars

SNPPath = "/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/Paper_kleinesPaper/DNASeq/";
SNPName = ["20210623_Wns11-19_SNPSummary", "20210623_Vns01-08_SNPSummary", "20210623_BAns01-08_SNPSummary", "20210623_Control_SNPSummary"];

allSNPSummaries = [];
for i = 1 : length(SNPName)
Data = load(SNPPath + SNPName{i} + ".mat");
        allSNPSummaries = [allSNPSummaries Data.SNPSummary];
end


CNPName = ["20210623_Wns11-19_CNPSummary", "20210623_Vns01-08_CNPSummary.mat", "20210623_BAns01-08_CNPSummary", "20210623_Control_CNPSummary"];

allCNPSummaries = [];

for i = 1 : length(CNPName)
Data = load(SNPPath + CNPName{i});
        allCNPSummaries = [allCNPSummaries Data.CNPSummary'];
end
