degQ = [3257091:3257600];

for i=1:numel(CNPSummary)
    clear infos
    if isempty(CNPSummary(i).Adist)
        continue
    end
    infos        = [CNPSummary(i).Adist(:,[2:3])];
    allPos_inClu = [];
    for cl = 1:size(infos,1)
        allPos_inClu = [allPos_inClu infos(cl,1):infos(cl,2)];
    end
    allPos_inClu_collect{i} = allPos_inClu;
end

% test for degQ

for i=1:numel(allPos_inClu_collect)
    overlap_idx{i} = find(ismember(degQ',allPos_inClu_collect{i}));
end