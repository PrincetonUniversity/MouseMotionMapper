function K = makePartitionPlot(L,partitionIdx,stateDefs,fillSpaces)

    
    if nargin < 3 || isempty(fillSpaces)
        fillSpaces = true;
    end


    q = unique(partitionIdx);
    N = length(q);
    s = size(L);
    
    useStateDefs = nargin > 2 && length(partitionIdx) == length(stateDefs);
    
    K = zeros(s);
    for i=1:N
        idx = partitionIdx == q(i);
        if useStateDefs
            K(ismember(L,stateDefs(idx))) = i;
        else
            K(ismember(L,idx)) = partitionIdx(i);
        end
    end
    
    
    if nargin > 3
        if fillSpaces
            [ii,jj] = find(L == 0);
            idx = ii > 1 & jj > 1 & ii < s(1) & jj < s(2);
            ii = ii(idx);
            jj = jj(idx);
            
            iAdds = [1;1;1;0;0;-1;-1;-1];
            jAdds = [1;0;-1;1;-1;1;0;-1];
            
            vals = zeros(size(ii));
            for i=1:length(vals)
                q = K(sub2ind(size(K),ii(i)+ iAdds,jj(i) + jAdds));
                vals(i) = mode(q(q>0));
            end
            
            K(sub2ind(size(K),ii,jj)) = vals;
            
        end
    end