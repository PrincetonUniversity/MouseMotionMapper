function [cGuesses,pClusts] = reembedHDK(projections, parameters, distEmbeddingData, C, K)

numModes = parameters.numModes;
[data,f] = findWavelets(projections,numModes,parameters);
pWave = log(data); pWave(data<-3)=-3;
% distEmbeddingData = log(distEmbeddingData);
% distEmbeddingData(distEmbeddingData<-3)=-3;

%pWave(:) = bsxfun(@rdivide,data,sum(data,2));
clear data

batchSize = 20000;
N = length(pWave(:,1));
batches = ceil(N/batchSize);
readout = 2000;
cGuesses = zeros(N,2);
pClusts = zeros(N,K);
for j = 1:batches
    fprintf(1,'\t Processing batch #%4i out of %4i\n',j,batches);
    idx = (1:batchSize) + (j-1)*batchSize;
    idx = idx(idx <= N);
    current_guesses = zeros(length(idx),2);
    currentData = pWave(idx,:);
	
    fprintf(1,['Size of current Data = ' num2str(size(currentData)) '\n']);

    D2 = pdist2(currentData,distEmbeddingData,'euclidean');
    
    [Dt Dt2] = sort(D2,2,'ascend');
    Dt = Dt(:,1:K); Dt2 = Dt2(:,1:K); DtC = C(Dt2);
    mddt = mode(DtC,2);
    current_guesses(:,1) = mddt;
    current_guesses(:,2) = mean(DtC==(repmat(mddt,1,K)),2);
    
    cGuesses(idx,:) = current_guesses;
    pClusts(idx,:) = DtC;
end

