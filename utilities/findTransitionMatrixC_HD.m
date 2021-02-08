function [bigallT,allTallDensities,bigAll_flux,allTs,allTdensities,allTstateDefs] = findTransitionMatrixC_HD(WR,eliminateZeros,eliminateSelfTransitions,M)

allStates = cell(size(WR));
for i = 1:length(WR)
    CC = bwconncomp(WR{i}==0);
    a = WR{i};
    for j = 1:CC.NumObjects
        if CC.PixelIdxList{j}(end)==length(a)
            a(CC.PixelIdxList{j}) = a(CC.PixelIdxList{j}(1)-1);
        else
            a(CC.PixelIdxList{j}) = a(CC.PixelIdxList{j}(end)+1);
        end
    end
    allStates{i} = a;
end

% M = max(max(LL2));
allTs = cell(size(allStates,1));
allTdensities = cell(size(allStates,1));
allTstateDefs = cell(size(allStates,1));
for i = 1:length(allStates)
    [allTs{i},allTdensities{i},allTstateDefs{i}]=findTransitionMatrix(allStates{i},eliminateZeros,eliminateSelfTransitions);
end
bigallT = zeros(M);
allTallDensities = zeros(M,1);
for i = 1:length(allStates)
    F = bsxfun(@times,allTs{i},allTdensities{i});
    bigallT(allTstateDefs{i},allTstateDefs{i}) = bigallT(allTstateDefs{i},allTstateDefs{i}) + F;
    allTallDensities(allTstateDefs{i}) = allTallDensities(allTstateDefs{i}) + allTdensities{i};
end

%bigT = bigT(allStateDefs,allStateDefs);

bigallT = bsxfun(@rdivide,bigallT,sum(bigallT,2));
allTallDensities = allTallDensities ./ sum(allTallDensities);
bigallT(isnan(bigallT)) = 0;
bigAll_flux = bsxfun(@times,bigallT,allTallDensities);



    