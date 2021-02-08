function [groups,segments,sizes] = makeGroupsAndSegments(watershedRegions,maxNum,startFrames,minLength)


    N = length(watershedRegions);

    segments = cell(maxNum,N);
    sizes = zeros(maxNum,N);
    for j=1:N
        for i=1:maxNum
            %CC = bwconncomp(watershedRegions{j}==i);
            CC = largeBWConnComp(watershedRegions{j} == i,minLength);
            sizes(i,j) = CC.NumObjects;
            segments{i,j} = zeros(CC.NumObjects,2);
            for k=1:CC.NumObjects
                segments{i,j}(k,:) = CC.PixelIdxList{k}([1 end]);
            end
        end
    end
    
    
    groups = cell(maxNum,1);
    for i=1:maxNum
        totalNum = sum(sizes(i,:));
        groups{i} = zeros(totalNum,3);
        count = 1;
        for j=1:N
            if ~isempty(segments{i,j})
                n = length(segments{i,j}(:,1));
                groups{i}(count + (0:(n-1)),:) = [zeros(n,1)+j segments{i,j}+startFrames(j)-1];
                count = count + n;
            end
        end
    end