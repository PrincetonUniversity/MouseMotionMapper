function p1Center = findP1Nose(p1)

tailP = double(squeeze(p1(16,:,:)));
centerP = double(squeeze(p1(1,:,:)));

tP = zeros(18,2,size(p1,3));
tP(:,1,:) = repmat(tailP(1,:),[18,1,1]);
tP(:,2,:) = repmat(tailP(2,:),[18,1,1]);
tP1 = double(p1)-double(tP);

tP2 = zeros(size(tP1));
for i  = 1:length(tailP)
    [tP2(:,1,i), tP2(:,2,i)] = cart2pol(tP1(:,1,i),tP1(:,2,i));
end

rotDir = zeros(length(tailP),1);
for i = 1:length(tailP)
    xtemp = tailP(:,i)-centerP(:,i);
    rotDir(i) = deg2rad(rem(atan2d(xtemp(1),xtemp(2))+360,360));
end

tP3 = zeros(size(tP2));
for i = 1:length(tailP)
    degtemp = tP2(:,1,i)+repmat(rotDir(i),[18 1]);
    [tP3(:,1,i), tP3(:,2,i)] = pol2cart(degtemp,tP2(:,2,i));
end

tP4 = zeros(size(tP3));
for i = 1:18
    tP4(i,1,:) = smooth(tP3(i,1,:),5);
    tP4(i,2,:) = smooth(tP3(i,2,:),5);
end

% tP4 = tP4(xIdx,:,:);
p1Center = reshape(tP4,[36 size(tP3,3)]);