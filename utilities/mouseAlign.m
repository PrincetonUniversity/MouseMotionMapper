function [] = mouseAlign(fileName,savePathName)

% dataName = '/Volumes/Mouse_Data/OFT_day2_0254/OFT-0311-00_box_data.mat';
% jointName = '/Volumes/Mouse_Data/OFT_day2_0254/OFT-0311-00_box_PREDICTED.h5';
[a b] = fileparts(fileName);

addpath(genpath('/home/klibaite/leap/leap/toolbox/'));
addpath(genpath('/home/klibaite/leap/leap/toolbox/utilities/'));
addpath(genpath('/home/klibaite/leap/leap/toolbox/aliases/'));
addpath(genpath('/home/klibaite/leap/leap/toolbox/hdf5/'));
addpath(genpath('/home/klibaite/leap/leap/toolbox/hdf5/hdf5prop/'));

dataName = [fileName(1:end-3) '_data.mat'];
jointName = [fileName(1:end-3) '_PREDICTED.h5'];

savePathInfo = [savePathName b '_aligned_info.mat'];
savePathAligned = [savePathName b '_aligned.h5'];

J = h5read(jointName,'/positions_pred');
load(dataName)
indexF = runDataBox.indexF;
m_s = runDataBox.s;
maxN = length(indexF);
framesN = runDataBox.framesN;
centroidsF = zeros(framesN,2);

for tt = 1:maxN
    fprintf(1,['Alignment round ' num2str(tt) ' of ' num2str(maxN) '\n']);
    currL = length(indexF(tt,1):indexF(tt,2));
    tailV = zeros(currL,2);
    rotDir = zeros(currL,1);
    tic
    I = h5readframes(fileName,'/box',indexF(tt,1):indexF(tt,2));
    toc
    I = squeeze(I(:,:,1,:));
    fprintf(1,['Read frames for round ' num2str(tt) ' of ' num2str(maxN) '\n']);
    % 1. Nose, 2. tailPt, 3. tail tip
    
    tempJoints = J(:,:,indexF(tt,1):indexF(tt,2));
    noseJ = double(squeeze(tempJoints(1,:,:)))';
    tailJ = double(squeeze(tempJoints(2,:,:)))';

    tempCent = mouseDataBox.centroid(indexF(tt,1):indexF(tt,2),:);
    
    noseCorr = noseJ+tempCent;
    tailCorr = tailJ+tempCent;

    noseCorrSmooth(:,1) = smooth(noseCorr(:,1),30); noseCorrSmooth(:,2) = smooth(noseCorr(:,2),50); 
    tailCorrSmooth(:,1) = smooth(tailCorr(:,1),30); tailCorrSmooth(:,2) = smooth(tailCorr(:,2),10);

    
    for i = 1:currL
        tailV(i,:) = [noseCorrSmooth(i,1)-tailCorrSmooth(i,1) (m_s(1)-noseCorrSmooth(i,2))-(m_s(1)-tailCorrSmooth(i,2))];
    end
    for i = 1:currL
        rotDir(i) = rem(atan2d(tailV(i,1),tailV(i,2))+360,360);
    end
    translation = tailCorrSmooth-tempCent-200;
    
    centroidsF(indexF(tt,1):indexF(tt,2),:) = tempCent+translation;
    translationF(indexF(tt,1):indexF(tt,2),:) = translation;
    rotVal(indexF(tt,1):indexF(tt,2),:) = rotDir;
    tic
    m_CX = zeros(400,400,currL);
    for i = 1:currL
        m_CX(:,:,i) = imrotate(imtranslate(I(:,:,i),[-translation(i,1) -translation(i,2)]),rotDir(i),'crop');
    end
    toc
    fprintf(1,['Generated images for round ' num2str(tt) ' of ' num2str(maxN) '\n']);
%     fprintf(1,['Writing Batch ' num2str(tt) '\n']);
%     flyMovie1 = VideoWriter(['flies_aligned_fly1_' num2str(tt) ]);
%     open(flyMovie1);
%     for f = 1:size(m_CX,3)
%         writeVideo(flyMovie1,uint8(m_CX(:,:,f)));
%     end
%     close(flyMovie1);
    tic
    outImages1 = reshape(m_CX,[400 400 1 currL]);
    clear m_CX
    if tt == 1
        dsetSize = [400 400 1 framesN];
        dataType = 'uint8';
        h5create(savePathAligned, '/box',dsetSize,'datatype',dataType);
        h5write(savePathAligned,'/box',uint8(outImages1),[1 1 1 indexF(tt,1)],[400 400 1 size(outImages1,4)])
    else
        h5write(savePathAligned,'/box',uint8(outImages1),[1 1 1 indexF(tt,1)],[400 400 1 size(outImages1,4)])
    end
    toc
    fprintf(1,['Wrote to file round ' num2str(tt) ' of ' num2str(maxN) '\n']);
    clear tempJoints noseJ tailJ noseCorrSmooth tailCorrSmooth noseCorr tailCorr tempCent tailV rotDir outImages1
end


mouseInfo.centroidsF = centroidsF;
mouseInfo.translationF = translationF;
mouseInfo.rotVal = rotVal;
runDataBox.jointPath = jointName;
runDataBox.savePathAligned = savePathAligned;
save(savePathInfo,'mouseInfo','parametersBox','runDataBox');
