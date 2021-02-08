function [] = runCluster(K,saveName,fileName)

addpath(genpath('./utilities/'));
addpath(genpath('./tSNE/'));

load(fileName);
CC = kmeans(trainingSetData,str2num(K),'Replicates',20,'Distance','cityblock','MaxIter',100);


save([saveName num2str(K) '_clusters.mat'],'K','CC','fileName');

