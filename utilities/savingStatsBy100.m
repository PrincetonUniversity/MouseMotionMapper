
parameters = setRunParameters([]);
parameters.pcaModes = 36;
parameters.samplingFreq = 80;
parameters.minF = .25;
parameters.maxF = 20;
numModes = 36;
numW = 25;
BIG_W = cell(12,5);
L_W = cell(12,5);
minLength = 20;
maxNum = 100; N = 100;


jv = 1:18;
jv2 = [jv jv+18];
% N = 53; maxNum = 53;
for w = 1:60
    fprintf(['Projections from mouse ' num2str(w) '\n']);
    for q = 1:5
        p1 = WT_joints_all{w,q};

        p1Center = findP1Center(p1);
        p1C = p1Center(jv2,:);
        for i = 1:36
        p1C(i,:) = medfilt1(smooth(p1C(i,:),5),5);
        end
        
        [W1 a] = findWavelets(p1C',36,parameters);
        
        WR = cC1CON{7}{w,q};
        G1 = makeGroupsAndSegments({WR},maxNum,1,minLength);
        % P1 = projections(:,1:13)+(projections(:,14:end));
        % P2 = sqrt((diff(projections(:,1:13))).^2+(diff(projections(:,14:end))).^2);
        % [W1 ~] = findWavelets(p1C',32,parameters);
        powerW = zeros(numModes,length(W1));
        for i = 1:numModes
            idxI = (numW*(i-1)+1):(numW*(i-1)+numW);
            powerW(i,:) = sum(W1(:,idxI),2);
        end
        BIdx = cell(1,N);
        SortedW1 = cell(size(BIdx));
        
        for i = 1:N
            gtemp = G1{i};
            btemp = [];
            SW1t = cell(size(gtemp,1),1);
            for j = 1:size(gtemp,1)
                btemp = [btemp gtemp(j,2):gtemp(j,3)];
                SW1t{j} = W1(gtemp(j,2):gtemp(j,3),:);
            end
            BIdx{i} = btemp;
            SortedW1{i} = SW1t;
        end
        
        allW = cell(size(BIdx));
        allP = cell(size(BIdx));
        for i = 1:N
            allW{i} = W1(BIdx{i},:);
            allP{i} = p1C(:,BIdx{i});
        end
        
        
        allWM = zeros(36,25,N);
        l_w = zeros(1,N);
        for i = 1:N
            tempwm = mean(allW{i},1);
            l_w(i) = size(allW{i},1);
            allWM(:,:,i) = reshape(tempwm,[25 36])';
        end
        
        L_W{w,q} = l_w;
        BIG_W{w,q} = allWM;
    end
end



W_all = zeros(36,25,N);
L_all = zeros(1,N);
for i = 1:N
    TW = zeros(36,25);
    TL = 0;
    for M = 1:60
        tw = BIG_W{M}(:,:,i);
        tl = L_W{M}(i);
        if sum(sum(isnan(tw)))==0
            TW = TW+tw;
            TL = TL+tl;
        end
    end
    W_all(:,:,i) = TW./TL;
    L_all(i) = TL;
end

wxall = W_all(1:18,:,:);
wyall = W_all(19:end,:,:);

Wxy = sqrt(wxall.^2+wyall.^2);




tp = zeros(1,100);
for i = 1:100
    tp(i) = sum(sum(Wxy(:,:,i)));
end



save('statsBy100.mat','Wxy','v100','std100','e100','estd100');










