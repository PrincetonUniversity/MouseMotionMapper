function [bvid,centroid,areas] = mouseSeg3(vid,framesize,varargin)

% MOUSESEG3 Segments out an individual mouse from video files. Intended to
% be used within makeMouseBox3.m (for open-field data).
%
%   [BVID,CENTROID,AREAS] = MOUSESEG3(VID,FRAMESIZE) finds the CENTROID
%   of the mouse in every frame of VID and outputs that along with a
%   bounding box video BVID and the approximate AREAS of the mouse object.
%       VID is a 3D matrix for the preprocessed video.
%       FRAMESIZE is a vector giving the dimensions of a single frame.
%       BVID is a video with 400x400 frames centered around CENTROID.
%       CENTROID is a matrix containing the coordinates of the mouse over
%       the course of VID.
%       AREAS is a vector containing the approximate area of the mouse,
%       based on a calculated threshold.
%
%   [BVID,CENTROID,AREAS] = MOUSESEG3(VID,FRAMESIZE,P) uses the parameters
%   specified in P to calculate AREAS.
%       P is a structure containing video parameters, assembled earlier in
%       makeMouseBox3.m.
%
%   Note that this function is currently optimized for OFT cup videos,
%   which are 1280x1024, with a single visible mouse, and certain imaging
%   conditions. Some things are hard-coded for convenience and may need to
%   be tweaked if this function is to be adapted for other purposes.
%
%   AREAS is carried over from mouseSeg2 (Ugne Klibaite) and is
%   questionably useful, so this code doesn't treat it that carefully
%   beyond what already was being used.
%
%   Author: Xiaoting Sun

% Error checking

if nargin < 2
    error('Not enough inputs.')
elseif nargin > 3
    error('Too many inputs.')
end

usethresh = true;

if nargin == 3
    if ~isstruct(varargin{3})
        warning('Optional third input must be a structure with field minPixelVal.')
    elseif ~isfield(varargin{3},'minPixelVal')
        warning('Optional third input must be a structure with field minPixelVal.')
    else
        p = varargin{3};
        usethresh = false;
    end
end

if numel(size(vid)) ~= 3
    error('vid must be a 3D matrix representing a video.')
end
if numel(size(framesize)) ~= 2 || numel(framesize) ~= 2
    error('framesize must be a vector with 2 elements.')
end

nframes = size(vid,3);
rfac = 8; % resize factor

% (I)nitialize mouse location

inframes = 5;
imaxframe = 20;
icent = NaN(inframes,2); % centroids
i = 1; % number of frames captured + 1
j = 1; % frame being considered

objMax = prod(framesize/rfac) / 10;

while i <= inframes && j <= imaxframe
    
    % resize and smooth image
    % note: faster to resize then Gaussian filter
    im = imgaussfilt(imresize(vid(:,:,j),1/rfac),5);
    
    % if nonzero frame, proceed to processing
    if sum(im(:)) > 0
        
        % get object info
        rp = regionprops(im>0,'Centroid','PixelList');
        
        % if only one normal-sized object, save that centroid
        if numel(rp) == 1 && length(rp.PixelList) < objMax
            icent(i,:) = rp.Centroid;
            i = i+1;
        else
            
            % get max brightness of each object
            objB = zeros(numel(rp),1);
            for k = 1:numel(rp)
                objIdx = sub2ind(size(im),rp(k).PixelList(:,2),rp(k).PixelList(:,1));
                objB(k) = max(im(objIdx));
            end
            
            % reset threshold to look at brightest object(s)
            rp = regionprops(im > max(objB)*0.8,'Centroid','PixelList');
            
            % if one object left, save that centroid
            if numel(rp) == 1
                icent(i,:) = rp.Centroid;
                i = i+1;
            else
            
                % look at brightest object(s)
                myobj = objB == max(objB);
                
                % if one is uniquely bright, save that centroid
                if sum(myobj) == 1
                    icent(i,:) = rp(myobj).Centroid;
                    i = i+1;
                end
            
            end
        
        end
        
    end
    
    j = j+1;
    
end

% take median centroid as initial value
prev = median(icent,1,'omitnan');

% default to center of image if inconclusive
if isnan(prev)
    prev = framesize / rfac / 2;
end

% Store initialized coordinates for later
iprev = prev;

% Process video

% initialize variables
mvid = zeros([framesize/rfac,nframes]); % mini video
idxZ = false(nframes,1); % zero frames (without any objects)
idxD = false(nframes,1); % frames where jump is too big
% separate in case of future tweaking
pts = NaN(nframes,2);

maxD = 30; % distance allowed b/w two frames

% run through v1
for i = 1:nframes
    
    % resize, process, segment, find objects
    mvid(:,:,i) = imresize(vid(:,:,i),1/rfac);
    im = imgaussfilt(mvid(:,:,i),5);
    seg = (im > mean([mean(im(:)) max(im(:))]));
    rp = regionprops(seg,'Centroid');
    
    if isempty(rp) % zero frame
        idxZ(i) = true;
    else
        objs = reshape([rp.Centroid],2,[])';
        dst = sqrt(sum((objs-prev).^2,2));
        pts(i,:) = objs(dst==min(dst),:);
        if min(dst) > maxD % big jump
            idxD(i) = true;
        else
            prev = pts(i,:);
        end
    end
    
end

% Check for mistakes and patch

idxCheck = find(idxZ+idxD); % all weird frames

disp('Uncertain frames in this chunk:')
disp(num2str(idxCheck))

thresh = mean(mvid(:)) + 4*std(mvid(:));

for i = 1:sum(idxCheck ~= 0)
    j = idxCheck(i);
    seg = mvid(:,:,j) > thresh;
    rp = regionprops(seg,'Area','Centroid');
    if isempty(rp) % zero frame again
        pts(j,:) = pts(j-1,:);
    else % a) biggest object; b) same as above (closest object)
%         objA = [rp.Area]';
        objs = reshape([rp.Centroid],2,[])';
        if j <= 1
            dst = sqrt(sum((objs-iprev).^2,2));
        else
            dst = sqrt(sum((objs-pts(j-1,:)).^2,2));
        end
        if min(dst) > maxD
            pts(j,:) = pts(j-1,:);
        else
            pts(j,:) = objs(dst==min(dst),:);
        end
    end
end

% Screen for "flickers" and interpolate (maybe in future editions?)

% Collect outputs

centroid = pts*rfac;
bvid = zeros(400,400,nframes); % bbox
areas = NaN(nframes,1);

if usethresh
    athresh = thresh;
else
    athresh = p.minPixelVal;
end

for i = 1:nframes
    x = (1-200:200)+round(centroid(i,1));
    y = (1-200:200)+round(centroid(i,2));
    bvid(:,:,i) = boxIn2(vid(:,:,i),x,y,framesize,[400 400]);
    areas(i) = sum(sum(bvid(:,:,i)>athresh));
end

end