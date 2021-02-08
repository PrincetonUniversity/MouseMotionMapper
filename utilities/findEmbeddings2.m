function [zValues,outputStatistics] = ...
	findEmbeddings2(projections,trainingData,trainingEmbedding,parameters)

fprintf(1,'Finding Wavelets\n');
numModes = parameters.pcaModes;

[data,f] = findWavelets(projections,numModes,parameters);
data = log(data); data(data<-3) = -3;
% trainingData = log(trainingData); trainingData(trainingData<-3) = -3;

% trainingData = log(trainingData); trainingData(trainingData<1e-6) = 1e-6;

fprintf(1,'Finding Embeddings\n');
    [zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
        findTDistributedProjections_fminPOS(data,trainingData,...
        trainingEmbedding,[],parameters);
    
    
    
    outputStatistics.zCosts = zCosts;
    outputStatistics.f = f;
    outputStatistics.numModes = numModes;
    outputStatistics.zGuesses = zGuesses;
    outputStatistics.inConvHull = inConvHull;
    outputStatistics.meanMax = meanMax;
    outputStatistics.exitFlags = exitFlags;




