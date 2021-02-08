function [] = makeHDK_proj(fileName,trainingC,trainingSet,savePath,saveS,pcaModes)

addpath(genpath('./utilities/'));
addpath(genpath('./tSNE/'));

load(fileName,'projections');
fN1 = fileName;
load(trainingC);
load(trainingSet);

fna = strsplit(fN1,'_');
fnb = fna{2};

KK = 11;

parameters = setRunParameters([]);
parameters.pcaModes = pcaModes;
parameters.samplingFreq = 80;
parameters.minF = .25;
parameters.maxF = 20;
parameters.numModes = parameters.pcaModes;

cGuesses = reembedHDK(projections,parameters,trainingSetData(1:4:end,:),CC(1:4:end),KK);

[a b] = fileparts(fN1);

save([savePath b saveS ],'cGuesses','fN1');
