function findSampleDataWithWavelets(file,parameters)


    if nargin < 2 
        parameters = makeParametersStructure_tsne([]);
    else
        parameters = makeParametersStructure_tsne(parameters);
    end
    

    skipPoints = parameters.skipPoints;
    dt = parameters.dt;
    omega0 = parameters.omega0;
    numModes = parameters.numModes;
    firstMode = parameters.firstMode;
    numPeriods = parameters.numPeriods;
    minT = parameters.minT;
    maxT = parameters.maxT;
    perplexity = parameters.samplePerplexity;
    relTol = parameters.sampleRelTol;
    
    [dir,name,~] = fileparts(file);
    idx = find(name=='_',1,'last');
    saveFile = [dir '/' name(1:idx-1) '_cluster.mat'];
    
    Ts = minT.*2.^((0:numPeriods-1).*log(maxT/minT)/(log(2)*(numPeriods-1)));
    f = fliplr(1./Ts);
    
    fprintf(1,'Loading Projections\n');
    load(file,'projections');
    projections = projections(:,firstMode + (1:numModes) - 1);
    
    
    N = length(projections(:,1));
    nF = length(f);
    d = nF*numModes;
    
    fprintf(1,'Calculating Wavelets\n');
    data = zeros(N,d);
    for j=1:numModes
        data(:,(1:nF) + (j-1)*nF) = ...
            fastWavelet_morlet_convolution_parallel(projections(:,j),f,omega0,dt)';
    end
    
    signalData = data(skipPoints:skipPoints:end,:);
    signalAmps = sum(signalData,2);
    signalData = bsxfun(@rdivide,signalData,signalAmps);
    
    clear data projections
    
    fprintf(1,'Running t-SNE\n');
    [D,~] = findKLDivergences(signalData);
    [yData,betas,~,costs] = tsne_d(D, 2, perplexity, relTol, parameters);
    

    save(saveFile,'yData','signalData','betas','costs','file','parameters');
    



