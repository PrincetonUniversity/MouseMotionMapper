function outPoints = drawCurvedArrow(x,y,a,plotLine)

    if nargin < 4 || isempty(plotLine)
        plotLine = false;
    end

    if iscolumn(x)
        x = x';
    end
    
    if iscolumn(y)
        y = y';
    end

    diff = y-x;
    mid = .5*(y+x);
    d = norm(diff);
    normVec = [-diff(2) diff(1)];
    normVec = normVec ./ d;
    
    bezPoints = [x; mid + a*d.*normVec; y];
    
    outPoints = bezier_(bezPoints);
    
    if plotLine
        plot(outPoints(:,1),outPoints(:,2),'b-')
    end