%% PIPELINE
% calculate distances from body part positions

xIdx = [1 2 5 6 7 8 12 13 14 15 16];
yIdx = [1 2 5 6 7 8 12 13 14 15 16];
% here we track the nose/snout, limbs, and tail point


[X Y] = meshgrid(xIdx,yIdx);
X = X(:); Y = Y(:);
IDX = find(X~=Y);


savePath = '/Volumes/fileset-shaevitz/klibaite/mouseREPO/MARCH_2020/distSub/';
for m = 1:60
    m
    for j = 1:5
        p1 = WT_joints_all{m,j};
        p1Dist = zeros(121,size(p1,3));
        for i = 1:size(p1Dist,1)
            p1Dist(i,:) = returnDist(squeeze(p1(X(i),:,:)),squeeze(p1(Y(i),:,:)));
        end
        
        p1Dsmooth = zeros(size(p1Dist));
        for i = 1:size(p1Dist,1)
            p1Dsmooth(i,:) = medfilt1(smooth(p1Dist(i,:),5),5);
        end
        p1Dist = p1Dsmooth(IDX,:);
        save([savePath 'WT_mouse_' num2str(m) '_day_' num2str(j) '_Distance.mat'], 'p1Dist');
    end
end


%% do online PCA to find modes

[s t] = unix('find `pwd` -name "*Distance.mat"');
t = strsplit(t,'\n');
tt = t(1:end-1);

batchSize = 20000;
num = length(tt);
firstBatch = true;
currentImage = 0;

for j = 1:num
    fprintf(1,'\t File #%5i out of %5i\n',j,num);
    load(tt{j});
    if firstBatch
        firstBatch = false;
        if size(p1Dist,2)<batchSize
            cBatchSize = size(p1Dist,2);
            X = p1Dist';
        else
            cBatchSize = batchSize;
            X = p1Dist(:,randperm(size(p1Dist,2),cBatchSize))';
        end
        
        currentImage = batchSize;
        muv = sum(X);
        C = cov(X).*batchSize + (muv'*muv)./cBatchSize;
    else
        
        if size(p1Dist,2)<batchSize
            cBatchSize = size(p1Dist,2);
            X = p1Dist';
        else
            cBatchSize = batchSize;
            X = p1Dist(:,randperm(size(p1Dist,2),cBatchSize))';
        end
        
        tempMu = sum(X);
        muv = muv+tempMu;
        C = C + cov(X).*cBatchSize + (tempMu'*tempMu)./cBatchSize;
        currentImage = currentImage + cBatchSize;
        fprintf(1,['Current Image = ' num2str(currentImage) '\n']);
    end
end

%save('/Volumes/fileset-shaevitz/klibaite/mouseREPO/MARCH_2020/vecsVals_subSelect.mat','C','muv','currentImage');

L = currentImage;
muv = muv./L;
C = C./L - muv'*muv;

[vecs,vals] = eig(C);
vals = flipud(diag(vals));
vecs = fliplr(vecs);

load('/Volumes/fileset-shaevitz/klibaite/mouseREPO/MARCH_2020/vecsVals_subSelect.mat','C','muv','currentImage','vecs','vals');



%p2Dist = bsxfun(@minus,p1Dist,muv');
%projections = p2Dist'*vecs(:,1:numProjections);


dataAll = [];

numProjections = 10;
numModes = 10; pcaModes = numModes;
minF = .25; maxF = 20;

parameters = setRunParameters([]);
parameters.trainingSetSize = 80;
parameters.pcaModes = pcaModes;
parameters.samplingFreq = 80;
parameters.minF = minF;
parameters.maxF = maxF;
numPerDataSet = parameters.trainingSetSize;
numPoints = 5000;

for t = 1:length(tt)
    t
    fileName = tt{t};
    load(fileName,'p1Dist');
    p2Dist = bsxfun(@minus,p1Dist,muv');
projections = p2Dist'*vecs(:,1:numProjections);

[data,~] = findWavelets(projections,numModes,parameters);
amps = sum(data,2);

N = length(projections(:,1));
numModes = parameters.pcaModes;
skipLength = floor(N / numPoints);
if skipLength == 0
    skipLength = 1;
    numPoints = N;
end

firstFrame = mod(N,numPoints) + 1;
signalIdx = firstFrame:skipLength:(firstFrame + (numPoints-1)*skipLength);
%signalData = bsxfun(@rdivide,data(signalIdx,:),amps(signalIdx));
signalAmps = amps(signalIdx);
nnData = log(data(signalIdx,:));
nnData(nnData<-3) = -3;

yData = tsne(nnData);
[signalData,signalAmps] = findTemplatesFromData(...
                    nnData,yData,signalAmps,numPerDataSet,parameters);

dataAll = [dataAll; signalData];

end

tic
ydata = tsne(dataAll);
toc
tic
C100 = kmeans(dataAll,100,'Replicates',10);
toc
tic
C200 = kmeans(dataAll,200,'Replicates',10);
toc

load('trainingSet_new10.mat','dataAll','ydata','C100','C200');




eV = embeddingValues;
eV(~outputStats.inConvHull,:) = outputStats.zGuesses(~outputStats.inConvHull,:);
scatter(eV(:,1),eV(:,2),'.'); axis([-60 60 -60 60])


fileName = '/Volumes/fileset-shaevitz/klibaite/mouseREPO/MARCH_2020'


addpath('/Volumes/Fly_Courtship/MultiFlyAnalysis/utilities/')
addpath('/Volumes/Fly_Courtship/MultiFlyAnalysis/tSNE/')
addpath('/Volumes/Fly_Courtship/MultiFlyAnalysis/visualization/')
addpath('/Volumes/Fly_Courtship/MultiFlyAnalysis/')


load(fileName,'p1Dist');
load(tData); %,'trainingSetData','ydata','cdata1','cdata2');
trainingSetData = dataAll; %trainingSetData(1:2:end,:);
trainingEmbedding = ydata; %ydata(1:2:end,:);

cdata1 = C100; % cdata1(1:2:end);
cdata2 = C200; % cdata2(1:2:end);

%trainingEmbedding = ydata;

fprintf(1,['training data = ' num2str(size(trainingSetData)) '\n']);
load('/Volumes/fileset-shaevitz/klibaite/mouseREPO/MARCH_2020/vecsVals_subSelect.mat');

numProjections = 10; numModes = 10; pcaModes = 10;

parameters = [];
parameters = setRunParameters([]);
parameters.pcaModes = 10;
parameters.samplingFreq = 80;
parameters.minF = .25;
parameters.maxF = 20;
parameters.numModes = numModes;

[s t] = fileparts(fileName);
fN = t;

p2Dist = bsxfun(@minus,p1Dist,muv');
projections = p2Dist'*vecs(:,1:numProjections);

fprintf(1,'Finding Wavelets\n');
numModes = parameters.pcaModes;

[data,f] = findWavelets(projections,numModes,parameters);
data = log(data); data(data<-3) = -3;
% trainingData = log(trainingData); trainingData(trainingData<-3) = -3;
trainingData = dataAll; trainingEmbedding = ydata;
% trainingData = log(trainingData); trainingData(trainingData<1e-6) = 1e-6;

fprintf(1,'Finding Embeddings\n');
    [zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
        findTDistributedProjections_fmin(data,trainingData,...
        trainingEmbedding,[],parameters);
    
    z = zValues; z(~inConvHull,:) = zGuesses(~inConvHull,:);

[embeddingValues,outputStats] = findEmbeddings2(projections,trainingSetData,trainingEmbedding,parameters);          
save(['/tigress/SHAEVITZ/klibaite/mouseREPO/MARCH_2020/EVnew/' fN '_RE.mat'],'embeddingValues','outputStats','fileName');

KK = 101;
[cGuesses,pClusts] = reembedHDK(projections,parameters,trainingSetData,cdata1,KK);
save(['/tigress/SHAEVITZ/klibaite/mouseREPO/MARCH_2020/C100new/' fN '_RE.mat'],'cGuesses','fileName','pClusts');
clear cGuesses pClusts

[cGuesses,pClusts] = reembedHDK(projections,parameters,trainingSetData,cdata2,KK);
save(['/tigress/SHAEVITZ/klibaite/mouseREPO/MARCH_2020/C200new/' fN '_RE.mat'],'cGuesses','fileName','pClusts');


%% read and organize

savePath = '/Volumes/fileset-shaevitz/klibaite/mouseREPO/MARCH_2020/';
s1 = 'EVnew/'; s2 = 'C100new/'; s3 = 'C200new/';

for i = 1:60
    for j = 1:5
        WTs1{i,j} = [savePath s1 'WT_mouse_' num2str(i) '_day_' num2str(j) '_Distance_RE.mat'];
        WTs2{i,j} = [savePath s2 'WT_mouse_' num2str(i) '_day_' num2str(j) '_Distance_RE.mat'];
        WTs3{i,j} = [savePath s3 'WT_mouse_' num2str(i) '_day_' num2str(j) '_Distance_RE.mat'];
    end
end

for i = 1:91
    for j = 1:4
        ASDs1{i,j} = [savePath s1 'ASD_' ASD_numbers{i,j} '_' ASD{i,j} '_' ASD_zygosity{i,j} '_day_' num2str(j) '_Distance_RE.mat'];
        ASDs2{i,j} = [savePath s2 'ASD_' ASD_numbers{i,j} '_' ASD{i,j} '_' ASD_zygosity{i,j} '_day_' num2str(j) '_Distance_RE.mat'];
        ASDs3{i,j} = [savePath s3 'ASD_' ASD_numbers{i,j} '_' ASD{i,j} '_' ASD_zygosity{i,j} '_day_' num2str(j) '_Distance_RE.mat'];
    end
end


    
load('/Volumes/fileset-shaevitz/klibaite/mouseREPO/MARCH_2020/e_filenames.mat','WTs1','WTs2','WTs3','ASDs1','ASDs2','ASDs3');

WTev = cell(size(WTs1)); WTC1 = WTev; WTC2 = WTev;
for i = 1:60
    fprintf(1,['Processing Mouse #' num2str(i) '\n']);
    for j = 1:5
        load(WTs1{i,j})
        eV = embeddingValues;
        eV(~outputStats.inConvHull,:) = outputStats.zGuesses(~outputStats.inConvHull,:);
        WTev{i,j} = eV;
        WTincv(i,j) = mean(outputStats.inConvHull);
        
        load(WTs2{i,j}); 
        WTC1{i,j} = cGuesses(:,1);
        
        load(WTs3{i,j});
        WTC2{i,j} = cGuesses(:,1);
    end
end
ASDev = cell(size(ASDs1)); ASDC1 = ASDev; ASDC2 = ASDev;
for i = 37:91
    fprintf(1,['Processing Mouse #' num2str(i) '\n']);
    for j = 1:4
        load(ASDs1{i,j})
        eV = embeddingValues;
        eV(~outputStats.inConvHull,:) = outputStats.zGuesses(~outputStats.inConvHull,:);
        ASDev{i,j} = eV;
        ASDincv(i,j) = mean(outputStats.inConvHull);
        
        load(ASDs2{i,j}); 
        ASDC1{i,j} = cGuesses(:,1);
        
        load(ASDs3{i,j});
        ASDC2{i,j} = cGuesses(:,1);
    end
end


save('tempSaveMarch12_allev.mat','WTev','WTC1','WTC2','ASDev','ASDC1','ASDC2');

for j = 1:5
    subplot(1,5,j);
tmp = WTev{1,j};
[xx,D] = findPointDensity(tmp,1.5,501,[-55 55 -55 55]);
imagesc(D); colormap(cmap1);
end


ev = [combineCells(WTev(:)); combineCells(ASDev(:))];
[xx, D] = findPointDensity(ev,.8,501,[-55 55 -55 55]);
LL = watershed(-D,8);
LL(D<=1e-6) = 0;


tmp = WTev{1,1};
tc = WTC2{1,1};
scatter(tmp(:,1),tmp(:,2),[],tc,'filled'); colormap(cmap1);

eV = embeddingValues;
eV(~outputStats.inConvHull,:) = outputStats.zGuesses(~outputStats.inConvHull,:);


load('/Volumes/Mouse_Data/colormapsMouse.mat')

eV = embeddingValues;
eV(~outputStats.inConvHull,:) = outputStats.zGuesses(~outputStats.inConvHull,:);

idx = 70000:80000;
plot(eV(idx,1),eV(idx,2),'k'); hold on
scatter(eV(idx,1),eV(idx,2),200,1:length(idx),'.')
colormap(cmap1)


% visualize by length of dwell?


load('/Volumes/Mouse_Data/mouse-behavior/ASD_WT_v.mat','vWT','vASD');


te = WTev{1,1};
tv = vWT{1,1};
scatter(te(:,1),te(:,2),[],tv,'filled'); axis ij; caxis([0 10])




fitOnly = 1;
numGMM = 2;
obj = [];
minRest = 5;
pThreshold = [];
medianLength = 2;
vSmooth = [];

[watershedRegions,segments,v,obj,pRest,vals,vx,vy] = ...
            findWatershedRegions_v2(ev,xx,LL,vSmooth,medianLength,pThreshold,minRest,obj,fitOnly,numGMM);
        
fitOnly = false;
for i = 1:60
    for j = 1:5
        [WTwr{i,j},segments,v,obj,pRest,vals,vx,vy] = ...
            findWatershedRegions_v2(WTev{i,j},xx,LL,vSmooth,medianLength,pThreshold,minRest,obj,fitOnly,numGMM);
    end
end

idx = 1:1000;
wr = WTwr{1,1};
figure(1); scatter(te(idx,1),te(idx,2),[],wr(idx),'filled');
colormap(cmap3)
figure(2); scatter(te(idx,1),te(idx,2),[],WTC1{1,1}(idx),'filled');
colormap(cmap3)














