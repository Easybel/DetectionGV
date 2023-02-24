InsertV = load('/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/Paper_kleinesPaper/DNASeq/20210625_Vns01-08_Insertions.mat'); 
InsertV = InsertV.insertions;
DelV = load('/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/Paper_kleinesPaper/DNASeq/20210625_Vns01-08_DelDup.mat'); 
DelV = DelV.deldup.del;

InsertW = load('/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/Paper_kleinesPaper/DNASeq/20210625_Wns11-19_Insertions.mat'); 
InsertW = InsertW.insertions;
DelW = load('/home/isabel/Dokumente/P1_EvolExp_CLASSIC_Bacillus/Paper_kleinesPaper/DNASeq/20210625_Wns11-19_DelDup.mat'); 
DelW = DelW.deldup.del;

maskInsertV10 = contains([InsertV.Sample],"10");
maskInsertW10 = contains([InsertW.Sample],"10");
maskDelV10 = contains([DelV.sample],"10");
maskDelW10 = contains([DelW.sample],"10");
InsertV10  = InsertV(maskInsertV10);
InsertW10  = InsertW(maskInsertW10);
DelV10 = DelV(maskDelV10);
DelW10 = DelW(maskDelW10);



maskInsertV20 = contains([InsertV.Sample],"20");
InsertV20  = InsertV(maskInsertV20);
maskInsertW20 = contains([InsertW.Sample],"20");
InsertW20  = InsertW(maskInsertW20);

figure(1)
p1 = plot([mean(transW10) mean(transV10) mean(transB10)],[sum([InsertW10.Length]) sum([InsertV10.Length]) 0],'*'); hold on;
p2 = plot([mean(transW20) mean(transV20) mean(transB20)],[sum([InsertW20.Length]) sum([InsertV20.Length]) 0],'*');