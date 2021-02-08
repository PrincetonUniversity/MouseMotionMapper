function out = combineCells(x,dim)

    x = x(returnCellLengths(x) > 0);

    L = length(x);
    s = size(x{1});
    if nargin < 2 || isempty(dim)
        [~,dim] = max(s);
    end
    
    if dim == 1
        
        d = s(2);
        lengths = returnCellLengths(x) ./ d;
        try
	cVals = [0; cumsum(lengths)];
        catch
	cVals = [0; cumsum(lengths)'];
	end
        out = zeros(cVals(end),d);
        
        for i=1:L
            out(cVals(i)+1:cVals(i+1),:) = x{i};            
        end
        
    else
        
        d = s(1);
        lengths = returnCellLengths(x) ./ d;
        cVals = [0 cumsum(lengths)];
        
        
        out = zeros(d,cVals(end));
        
        for i=1:L
            out(:,cVals(i)+1:cVals(i+1)) = x{i};            
        end
        
    end
