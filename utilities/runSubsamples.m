function [] = runSubsamples(fileName,pcaModes,minF,maxF,savePath)

addpath(genpath('./utilities/'));
addpath(genpath('./tSNE/'));

[s t] = fileparts(fileName);
fN = t;

parameters =  setRunParameters([]);
parameters.trainingSetSize = 400;
parameters.pcaModes = pcaModes;
parameters.samplingFreq = 80;
parameters.minF = minF;
parameters.maxF = maxF;
numPerDataSet = parameters.trainingSetSize;

[yData,signalData,signalAmps,~] = file_embeddingSubSampling(fileName,parameters);    
[trainingSetData,trainingSetAmps] = findTemplatesFromData(signalData,yData,signalAmps,numPerDataSet,parameters);
                                                                                    
save([savePath t '_subsample.mat'],'trainingSetAmps','trainingSetData','fileName');

