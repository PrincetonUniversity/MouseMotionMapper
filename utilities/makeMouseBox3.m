function makeMouseBox3(filename,savePath)

% MAKEMOUSEBOX3 Segments out mouse from HDF5 video files.
%
%   MAKEMOUSEBOX3(FILEPATH,SAVEPATH) processes video specified by FILEPATH
%   and saves the results to SAVEPATH. Outputs are bounding box video with
%   400x400 frames, centroid coordinates, and approximate object areas.
%
%   Optimized for OFT videos, which are 1280x1024 and have certain imaging
%   conditions.
%
%   Essentially and edited/streamlined version of makeMouseBox2.
%
%   Author: Ugne Klibaite
%   Adapted by Xiaoting Sun

% Error checking

if ~ischar(filename) || ~ischar(savePath)
    error('Inputs must be strings')
end

% File details

[~,name,ext] = fileparts(filename);

if ~strcmp(ext,'.h5')
    error('Must be an HDF5 file.')
end

if ~strcmp(savePath(end),'/')
    savePath = [savePath '/'];
end

myfile = [savePath name '_box.h5'];
mydata = [savePath name '_box_data.mat'];

p.savePath = myfile;
p.savePathData = mydata;
p.fileName = filename;

fprintf(1,['Input: ' filename '\n']);
fprintf(1,['Output: ' myfile '\n']);

% Add paths

addpath(genpath('leap'));
addpath(genpath('leap/leap'));
addpath(genpath('trackAndAlign'));
addpath(genpath('trackAndAlign/utilities'));

% Set parameter

p.chunkSize = 10000;
p.minPixelVal = 5;
p.minTailVal = 30;

% Define run indices

h = h5info(filename,'/pg0');
p.dataName = h.Name;
nframes = h.Dataspace.Size(4);
maxN = ceil(nframes/p.chunkSize);

% start and end indices of each chunk
cIdx = [((1:maxN)-1)*p.chunkSize+1; (1:maxN)*p.chunkSize]';
if cIdx(maxN,2) > nframes
    cIdx(maxN,2) = nframes;
end

% Make median image

tic

% note: h5readframes isn't very efficient with big datasets + rand
% sampIdx = randsample(1:nframes,50);
sampIdx = randsample(1:round(nframes/50),1):round(nframes/50):nframes;
sampFrames = h5readframes(filename,p.dataName,sampIdx);

medIm = median(squeeze(sampFrames(:,:,1,:)),3);
s = size(medIm);
clear sampFrames

toc
fprintf(1,'Generated median image.\n');

% Find centroids

pts = zeros(nframes,2);
myareas = NaN(1,nframes);

for n = 1:maxN
    
    fprintf(1,['Segmentation round ' num2str(n) ' of ' num2str(maxN) '.\n']);
    
    tic
    bigMovie = h5read(filename,['/' p.dataName],[1 1 1 cIdx(n,1)],...
        [Inf Inf Inf cIdx(n,2)-cIdx(n,1)+1]);
    toc
    
    bigMovie = squeeze(bigMovie(:,:,1,:));
    z = size(bigMovie,3);
    
    tic
    [bbox,pts(cIdx(n,1):cIdx(n,2),:),myareas(cIdx(n,1):cIdx(n,2))] =  ...
    mouseSeg3(bigMovie-repmat(medIm,[1 1 size(bigMovie,3)]),s);
    toc
   
    bbox2 = reshape(bbox,[400 400 1 z]);
 
    if n == 1
        dsetSize = [400 400 1 nframes];
        dataType = 'uint8';
        h5create(p.savePath, '/box',dsetSize,'datatype',dataType);
        h5write([p.savePath],'/box',uint8(bbox2),[1 1 1 cIdx(n,1)],[400 400 1 size(bbox2,4)])
        
    else
        h5write([p.savePath],'/box',uint8(bbox2),[1 1 1 cIdx(n,1)],[400 400 1 size(bbox2,4)])
    end
    
    clear bbox bbox2 bigMovie

end

% Collect outputs

mouseDataBox.centroid = pts;
mouseDataBox.areas = myareas;

runDataBox.medIm = medIm;
runDataBox.fileName = filename;
runDataBox.savePath = myfile;
runDataBox.h = h;
runDataBox.framesN = nframes;
runDataBox.indexF = cIdx;
runDataBox.s = s;

parametersBox = p;

save(mydata,'mouseDataBox','runDataBox','parametersBox');

end