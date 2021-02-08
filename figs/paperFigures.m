%% Paper Figures - [Klibaite et al]

% colorbars to recreate figures
% Cntnap2(+/+), Cntnap2(+/-), Cntnap2(-/-), L7-Tsc1 NEG, L7-Tsc1 HET,
% L7-Tsc1 HOM, Black6 MALE, Black6 FEMALE
cvcon = [188 128 189; 103 169 207; 2 129 138; ...
    253 141 60; 227 26 28; 128 0 38; 0 0 0; 248 129 191]./256;


% idle, groom, slow exp, fast exp, rear, climb, turning, locomotion
cbarbeh = [26 26 255;
    102 205 255;
    225 179 217; %204 204 255;
    172 0 230; %152 78 163;
    102 255 153;
    0 204 0;
    255 204 0;
    228 26 28]./256;


% Notes for plotting:
% con always refers to condition as ordered above
% m refers to mouse number
% d refers to day
% b refers to behavior (often out of 8 classes)

%% Preprocessing
% behavioral labels applied to each timepoint (at 80 Hz) are stored in cell
% arrays 

% C1CON - cell array with 8 conditions, mouse x day for each condition -
% 100 clusters unsorted

% sC1CON - cell array with 8 conditions, mouse x day for each condition - 8
% classes sorted after manual curation

cd klibaite_et_al/

load('data/mouseSkeletonSimple_400.mat')
load('data/allC100_CON.mat')
load('data/allC100_8C_CON.mat')
load('data/jointsCON.mat');

% fill any NaN values
cC1CON = cell(size(sC1CON));
for con = 1:8
    [mN dN] = size(sC1CON{con});
    for m = 1:mN
        for d = 1:dN
            temp = sC1CON{con}{m,d};
            temp2 = fillmissing(temp,'nearest');
            cC1CON{con}{m,d} = temp2;
        end
    end
end

% fill any NaN values
fC1CON = cell(size(C1CON));
for con = 1:8
    [mN dN] = size(C1CON{con});
    for m = 1:mN
        for d = 1:dN
            temp = C1CON{con}{m,d};
            temp2 = fillmissing(temp,'nearest');
            fC1CON{con}{m,d} = temp2;
        end
    end
end

% derive CON_P (center point), and CON_V (speed)
for con = 1:8
    fprintf(['Projections from CON ' num2str(con) '\n']);
    jCON = jointsCON{con};
    [mN dN] = size(jCON);
    %cCON = centroidsCON{con};
    for m = 1:mN
        fprintf(['Projections from mouse ' num2str(m) '\n']);
        for d = 1:dN
            try
            ja20 = jCON{m,d}; n = size(ja20,3);
            
            nJoints = size(ja20,1);
            i = 9; jt1 = squeeze(ja20(i,:,:)); jt1(1,:) = smooth(medfilt1(jt1(1,:),4),8); jt1(2,:) = smooth(medfilt1(jt1(2,:),4),8);
            i = 16; jt2 = squeeze(ja20(i,:,:)); jt2(1,:) = smooth(medfilt1(jt2(1,:),4),8); jt2(2,:) = smooth(medfilt1(jt2(2,:),4),8);
            centpt = zeros(n,2); heada = zeros(n,1);
            centpt(:,1) = (jt1(1,:)+jt2(1,:))./2; centpt(:,2) = (jt1(2,:)+jt2(2,:))./2;
            dir = jt2-jt1;
            for i = 1:n
                heada(i) = atan2d(dir(2,i),dir(1,i));
            end
            vt = [0 0; diff(centpt)]; vt = sqrt(sum(vt.^2,2)); cent = vt;
            p1Center = findP1Center(ja20);
            
            CON_P{con}{m,d} = centpt;
            CON_V{con}{m,d} = vt;
            catch
            end
        end
    end
end

% or just load it here
%load('data/centroids_velocities_CON.mat')


%% figure 1
% sample timepoints used in figure 1 
% plot centroid trajectory

% joints and index to use - WT mouse 3 on day 1
juse = [1 8 7 12 13 16 17 18];
idx = 40001:49600;
wtemp = sC1CON{7}{3,1};
%wtempC = sortid8(wtemp);
wtemp = fillmissing(wtemp,'nearest');
centroids = CON_P{7}{3,1};
v = smooth(CON_V{7}{3,1},8)*(80/(1000*1.97));
jj = jointsCON{7}{3,1};

jointst = jj(:,:,idx);
centt = centroids(idx,:);
%ht = H(idx);
wt2 = wtemp(idx);


% plot joints in arena for timepoints - specify joints and colors
skel.joints_idx = [1 2; 2 9; 9 16; 16 17; 17 18; 1 3; 1 4; 5 7; 6 8; 12 14; 13 15];
skel.color = [0 .447 .741;
    0.85 0.32 .098;
    .929 0.694 .125;
    .466 .674 .188;
    .737 .502 .7412; %
    .635 .078 .184;
    .85 .325 .098;
     0 .447 .741;
    .494 .1840 .556;
    .466 .674 .188
    .301 .745 .933;];

dd = 1:18; dd([9 10 11]) = [];

hold on
plot(centt(:,1),centt(:,2),'Color','k','LineStyle','-'); hold on;
for i = [1 200 300 500 1020 1200]
    scatter(jointst(dd,1,i),jointst(dd,2,i),[],'k','filled'); hold on;
    for j = 1:11
        xx = jointst(skel.joints_idx(j,:),1,i)'; 
        yy = jointst(skel.joints_idx(j,:),2,i)';
        plot(xx,yy,'Color',skel.color(j,:),'LineWidth',3);
    end
end
axis equal off
% saveas(gcf,['/Users/ugne/Dropbox/julyFigs/posturesIDX.pdf'])


% plot ethogram
strand = wtemp(idx);
for i = 1:8
n = length(strand);
isl = find(strand==i); ww = zeros(size(strand)); ww(isl)=1;
cc = bwconncomp(ww);
plot([0 n],[i i],'Color','k','LineWidth',.5); hold on;
for j = 1:length(cc.PixelIdxList)
    y = [i i]; x = [cc.PixelIdxList{j}(1) cc.PixelIdxList{j}(end)];
    plot(x,y,'Color',cbarbeh(i,:),'LineWidth',10);
end
end
axis off
saveas(gcf,['/Users/ugne/Dropbox/julyFigs/fullotherethogramIDX.tiff'])



% plot joint wavelet sample
nosedir = squeeze(jj(1,:,:))';
nosedir(:,1) = medfilt1(nosedir(:,1),20); nosedir(:,2) = medfilt1(nosedir(:,2),20);

t1dir = squeeze(jj(17,:,:))';
t1dir(:,1) = medfilt1(t1dir(:,1),20); t1dir(:,2) = medfilt1(t1dir(:,2),20);
t2dir = squeeze(jj(18,:,:))';
t2dir(:,1) = medfilt1(t2dir(:,1),20); t2dir(:,2) = medfilt1(t2dir(:,2),20);

waves = cell(18,1);
numModes = 1; pcaModes = numModes;
samplingFreq = 80;
minF = .25; maxF = 20;
omega0 = 5; numPeriods = 25; dt = 1./80;
minT = 1./maxF; maxT = 1./minF;
Ts = minT.*2.^((0:numPeriods-1).*log(maxT/minT)/(log(2)*(numPeriods-1)));
f = fliplr(1./Ts);
dt = 1 ./ samplingFreq;
minT = 1 ./ maxF;
maxT = 1 ./ minF;
%N = n;
for i = 1:18
    jt = squeeze(jj(i,:,:)); jt1 = medfilt1(jt(1,:),8); jt2 = medfilt1(jt(2,:),8);
    amplitudes1 = fastWavelet_morlet_convolution_parallel(jt1(1,:)',f,omega0,dt);
    amplitudes2 = fastWavelet_morlet_convolution_parallel(jt2(1,:)',f,omega0,dt);
    waves{i} = (amplitudes1+amplitudes2)./2;
end

% waves for juse at idx - 2 minutes
juse2 = juse([1:5 8]);
widx = [];
for i = 1:6
    jid = juse2(i);
    tw = waves{jid}(:,idx);
    widx = [widx; tw];
end

imagesc(widx);
axis off; colorbar;
%saveas(gcf,['/Users/ugne/Dropbox/julyFigs/waveletsIDX.tiff'])

% plot velocity
plot(v(idx),'Color','k','LineWidth',1); axis off
%saveas(gcf,['/Users/ugne/Dropbox/julyFigs/vIDX.tiff'])


% plot joint trajectories in x or y aligned space
jrot = findP1Nose(jointst);
jrot2 = findP1Center(jointst);

jr = findP1Nose(jj);
jr2 = findP1Center(jj);

jrotx = jrot(1:18,:);
jrotx2 = jrot2(1:18,:);
jroty = jrot(19:end,:);
jrot2y = jrot2(19:end,:);

cj = [0 1 0; 0 0 1; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 255/255 153/255 0/255; 166/255 77/255 255/255];
for i = 1:8
    plot(smooth(jrotx2(juse(i),4401:5000),4),'LineWidth',2,'Color',cj(i,:)); hold on
end
axis ij off
%saveas(gcf,['/Users/ugne/Dropbox/julyFigs/groom_x.pdf'])


for i = 1:8
    plot(smooth(jroty(juse(i),901:1500),4),'LineWidth',2,'Color',cj(i,:)); hold on
end
axis ij off
%saveas(gcf,['/Users/ugne/Dropbox/julyFigs/loc_y.pdf'])



% power spectrum from joints - WT only example
parameters = setRunParameters([]);
parameters.pcaModes = 36;
parameters.samplingFreq = 80;
parameters.minF = .25;
parameters.maxF = 20;
numModes = 36;
numW = 25;
BIG_W = cell(12,5);
L_W = cell(12,5);
minLength = 20;
maxNum = 100; N = 100;

WT_joints_all = jointsCON{7}; WTC1 = C1CON{7};
jv = 1:18;
jv2 = [jv jv+18];
for w = 1:5:60
    fprintf(['Projections from mouse ' num2str(w) '\n']);
    for q = 1:5
        p1 = WT_joints_all{w,q};

        p1Center = findP1Center(p1);
        p1C = p1Center(jv2,:);
        for i = 1:24
        p1C(i,:) = medfilt1(smooth(p1C(i,:),5),5);
        end
        
        [W1 a] = findWavelets(p1C',36,parameters);
        
        WR = WTC1{w,q};
        G1 = makeGroupsAndSegments({WR},maxNum,1,minLength);
        % P1 = projections(:,1:13)+(projections(:,14:end));
        % P2 = sqrt((diff(projections(:,1:13))).^2+(diff(projections(:,14:end))).^2);
        % [W1 ~] = findWavelets(p1C',32,parameters);
        powerW = zeros(numModes,length(W1));
        for i = 1:numModes
            idxI = (numW*(i-1)+1):(numW*(i-1)+numW);
            powerW(i,:) = sum(W1(:,idxI),2);
        end
        BIdx = cell(1,N);
        SortedW1 = cell(size(BIdx));
        
        for i = 1:N
            gtemp = G1{i};
            btemp = [];
            SW1t = cell(size(gtemp,1),1);
            for j = 1:size(gtemp,1)
                btemp = [btemp gtemp(j,2):gtemp(j,3)];
                SW1t{j} = W1(gtemp(j,2):gtemp(j,3),:);
            end
            BIdx{i} = btemp;
            SortedW1{i} = SW1t;
        end
        
        allW = cell(size(BIdx));
        allP = cell(size(BIdx));
        for i = 1:N
            allW{i} = W1(BIdx{i},:);
            allP{i} = p1C(:,BIdx{i});
        end
        
        allWM = zeros(36,25,N);
        l_w = zeros(1,N);
        for i = 1:N
            tempwm = mean(allW{i},1);
            l_w(i) = size(allW{i},1);
            allWM(:,:,i) = reshape(tempwm,[25 36])';
        end
        
        L_W{w,q} = l_w;
        BIG_W{w,q} = allWM;
    end
end

W_all = zeros(36,25,N);
L_all = zeros(1,N);
for i = 1:N
    TW = zeros(36,25);
    TL = 0;
    for M = 1:5:60
        tw = BIG_W{M}(:,:,i);
        tl = L_W{M}(i);
        if sum(sum(isnan(tw)))==0
            TW = TW+tw;
            TL = TL+tl;
        end
    end
    W_all(:,:,i) = TW./TL;
    L_all(i) = TL;
end

wxall = W_all(1:18,:,:);
wyall = W_all(19:end,:,:);

Wxy = sqrt(wxall.^2+wyall.^2);
b = 2; % example, one cluster from fast locomotion class
juse1 = [1 8 7 12 13 17 18]; % snout, linbs, tail
imagesc(Wxy(juse1,:,b)./sum(sum(Wxy(juse1,:,b)))); axis off equal


for b = 1:100
    figure(b); 
    imagesc(Wxy(juse1,:,b)./sum(sum(Wxy(juse1,:,b)))); colormap(cmap1);
    caxis([0 .025]); axis off equal
    saveas(gcf,['MARCH_2020/figs_MARCH19/powers/pow_b_' num2str(b) '.tiff']);
    close(figure(b));
end

% spatial labeling - from figure 3 extras
hold on
scatter(centt(:,1),centt(:,2),[],cvcon(wtemp(idx),:),'filled'); hold on;



%% figure 2

% heatmap of behavioral class usage:
load('data/info_8con_8b.mat')
load('data/statsBy100.mat')
for b = 1:8
    fib{b} = find(sortid8==b);
end
sfib = cell(size(fib));
for bg = 1:8
    fix = fib{bg};
    [~,sfix] = sort(v100(fix));
    sid = fix(sfix);
    sfib{bg} = sid;
end
vsorted = [];
for i = 1:8
    t = sfib{i};
    vsorted = [vsorted v100(t)];
end
combS = [];
for i = 1:8
    combS = [combS sfib{i}];
end

bf = cell(1,8);
for con = 1:8
    ct = fC1CON{con};
    [mN dN] = size(ct);
    bmd = zeros(mN,100,dN);
    for m = 1:mN
        for d = 1:dN
            temp = ct{m,d};
            htemp = hist(temp,1:100);
            htf = htemp./(sum(htemp));
            bmd(m,:,d) = htf(combS);
        end
    end
    bf{con} = bmd;
end

con = 7; d = 1;
temp = squeeze(bf{con}(:,:,d));
imagesc(temp); 


% 8-class
bf = cell(1,8);
for con = 1:8
    ct = sC1CON{con};
    [mN dN] = size(ct);
    bmd = zeros(mN,8,dN);
    for m = 1:mN
        for d = 1:dN
            temp = ct{m,d};
            htemp = hist(temp,1:8);
            htf = htemp./(sum(htemp));
            bmd(m,:,d) = htf;
        end
    end
    bf{con} = bmd;
end

d = 1; con = 7;
imagesc(bf{con}(:,:,d)); 



%% figure 3

% basic measures - 
% time spent in corners, time spent in center, 

mask = 2*ones(1200,1200);
mask(1:420,1:300) = 3; mask(850:1200,1:300) = 3;
mask(1:420,750:1200) = 3; mask(850:1200,750:1200) = 3;
mask(420:850,300:750) = 1;


m1 = mask==1;
m2 = mask==2;
m3 = mask==3;

% check where points fall relative to corners/edges
con = 1;
figure(con);
[mN dN] = size(CON_P{con});
if mN>10
    mN = 10;
end
dN = 4;
cpwt = CON_P{con};
for m = 1:mN
    for d = 1:dN
        subplot(mN,dN,4*(m-1)+d);
        imagesc(mask); hold on; scatter(cpwt{m,d}(1:100:end,1),cpwt{m,d}(1:100:end,2),1,'w','.')
    end
end

tz = nan(8,4,3); % condition,day, zone
for con = 1:8
    for d = 1:4
        temp = combineCells(CON_P{con}(:,d));
        tempz = zeros(length(temp),1);
            
            for ii = 1:length(temp)
                tempz(ii) = mask(ceil(temp(ii,1)),ceil(temp(ii,2)));
            end
        tz(con,d,:) = hist(tempz,1:3)./length(tempz);
    end
end

subplot(1,3,1); imagesc(tz(:,:,1));
subplot(1,3,2); imagesc(tz(:,:,2));
subplot(1,3,3); imagesc(tz(:,:,3));

% change between day1 and 4 for condition
con = 6; zone = 3; 
(tz(con,4,3)-tz(con,1,3))./tz(con,4,3)
con = 5; zone = 3; 
(tz(con,4,3)-tz(con,1,3))./tz(con,4,3)
con = 4; zone = 3; 
(tz(con,4,3)-tz(con,1,3))./tz(con,4,3)


dataf = []; datafN = [];
for con = 1:8
    [mN dN] = size(CON_P{con});
    for m = 1:mN
        for d = 1:4
            tempc = CON_P{con}{m,d};
            tempb = cC1CON{con}{m,d};
            tempz = zeros(size(tempb));
            
            for ii = 1:length(tempb)
                tempz(ii) = mask(ceil(tempc(ii,1)),ceil(tempc(ii,2)));
            end
            
            bzmat = nan(8,3);
            for bi = 1:8
                for zi = 1:3
                    bzmat(bi,zi) = sum(tempz==zi & tempb==bi)./length(tempz);
                end
            end
            for b = 1:8
                datatemp = [con d b bzmat(b,:)];
                datatempN = [con d b bzmat(b,:)./sum(bzmat(b,:))];
                dataf = [dataf; datatemp];
                datafN = [datafN; datatempN];
            end
        end
    end
end

Cz = cell(8,8,4,3);
for b = 1:8
user = datafN(find(datafN(:,3)==b),:);
for con = 1:8
    Ct = user(user(:,1)==con,:);
    for d = 1:4
        Ctd = Ct(Ct(:,2)==d,:);
        for zone = 1:3
            Cz{con,b,d,zone} = Ctd(:,zone+3);
        end
    end
end
end

mask2 = mask;
mask2(mask==3) = 2;
crossings = cell(8,4);

for con = 1:8
    [mN dN] = size(CON_P{con});
    for d = 1:4
        clear cx
        for m = 1:mN
            tempc = CON_P{con}{m,d};
            tempb = cC1CON{con}{m,d};
            tempz = zeros(size(tempb));
            for ii = 1:length(tempb)
                tempz(ii) = mask2(ceil(tempc(ii,1)),ceil(tempc(ii,2)));
            end
            cx(m) = sum(diff(tempz)~=0);
        end
        crossings{con,d} = cx;
    end
end
crossings{8,2}(16) = NaN;

conO = [7 8 1 2 3 4 5 6];
C = cell(8,4); cv = []; mids = []; cv2 = [];
for con = 1:8
    cont = conO(con);
    mo = (con-1)*5;
    for d = 1:4
        C{con,d} = crossings{cont,d}';
        mids = [mids mo+d];
        cv = [cv; 0 0 0];
        cv2 = [cv2; cvcon(cont,:)];
    end
end
CD = C'; CD = CD(:); 
oi = myBoxNA(CD,mids,0.4,cv,addPoints);
axis([0 41 0 500]); hold on;
oim = oi(:,1);
for con = 1:8
    cont = conO(con);
    mo = (con-1)*5;
    ds = (1:4)+mo; cid = (1:4)+(con-1)*4;
    scatter(ds,oim(cid),[],cvcon(cont,:),'filled');
    plot(ds,oim(cid),'Color',cvcon(cont,:),'LineWidth',2);
end
set(gcf,'Units','Normalized','OuterPosition',[0 0.04 .8 0.5])
%saveas

% distace covered - supplemental info
for con = 1:8
    temp = CON_V{con};
    dcc = zeros(size(temp));
    [mN dN] = size(temp);
    for m = 1:mN
        for d = 1:dN
            cd1 = temp{m,d}; cd2 = cd1; cd2(cd1>10) = 10; cd2 = smooth(cd2,20);
            dcc(m,d) = sum(sum(cd2))./(1.92*1000); % meters covered
        end
    end
    distCON{con} = dcc;
end

for con = 1:8
    figure(con);
    dN = size(distCON{con},2);
    C = cell(1,dN);
    for d = 1:dN
        C{d} = distCON{con}(:,d);
    end
    mids = [1:dN]; w = .4; cv = repmat(cvcon(con,:),[dN 1]); addPoints = 'true';
    myBoxNA(C,mids,w,cv,addPoints);
    axis([0 6 0 160]);
    %saveas
    axis off
end

% median locomotion speed - supp
medvCON = cell(8,1);
for con = 1:8
    temp = CON_V{con};
    dcc = zeros(size(temp));
    [mN dN] = size(temp);
    for m = 1:mN
        for d = 1:dN
            tempb = cC1CON{con}{m,d}; tbx = find(tempb==8);
            cd1 = temp{m,d}; cd2 = cd1; cd2(cd1>10) = 10; cd2 = smooth(cd2,20);
            dcc(m,d) = mean(cd2(tbx))*80/(1.92*1000); % meters covered
        end
    end
    medvCON{con} = dcc;
end

conO = [7 8 1 2 3 4 5 6];
C = cell(8,4); cv = []; mids = []; cv2 = [];
for con = 1:8
    cont = conO(con);
    mo = (con-1)*5;
    for d = 1:4
        C{con,d} = medvCON{cont}(:,d)';
        mids = [mids mo+d];
        cv = [cv; 0 0 0];
        cv2 = [cv2; cvcon(cont,:)];
    end
end

CD = C'; CD = CD(:); 
oi = myBoxNA(CD',mids,w,cv,addPoints);
axis([0 41 0 .5]); hold on;
oim = oi(:,1);
for con = 1:8
    cont = conO(con);
    mo = (con-1)*5;
    ds = (1:4)+mo; cid = (1:4)+(con-1)*4;
    scatter(ds,oim(cid),[],cvcon(cont,:),'filled');
    plot(ds,oim(cid),'Color',cvcon(cont,:),'LineWidth',2);
end
set(gcf,'Units','Normalized','OuterPosition',[0 0.04 .8 0.5])


(mean(C{7,1})-mean(C{8,1}))./mean(C{7,1})

(mean(C{6,1})-mean(C{5,1}))./mean(C{5,1})
(mean(C{6,1})-mean(C{4,1}))./mean(C{4,1})

%% figure 4 + 5
% 95% confidence interval for behavioral class
% choose window to average over, and how often to sample
% behavioral class by day
% sliding window over specific time-bin

w2 = 95998; L = 8;
window = 120*80;
sampleR = 10*80;
starts = 1:sampleR:(w2-window);
mids = (window/2):sampleR:(w2-(window/2));
ends = window:sampleR:w2;
BBAll = cell(8,1);

for con = 1:8
    con
    sortedTemp = cC1CON{con};
    mouseN = size(sortedTemp,1);
    dayN = size(sortedTemp,2);
    w2 = length(sortedTemp{1,1});
    BB = cell(L,dayN);
    for d = 1:dayN
        hh = zeros(L,length(starts),mouseN);
        for t = 1:length(starts)
            
            tempStatesIdx = starts(t):ends(t);
            for i = 1:mouseN
                try
                tempStates = sortedTemp{i,d}(tempStatesIdx);
                hh(:,t,i) = hist(tempStates,1:L)./length(tempStates);
                catch
                    hh(:,t,i) = nan(L,1);
                end
            end
        end
        for i = 1:L
            BB{i,d} = squeeze(hh(i,:,:))';
        end
    end
    BBAll{con,1} = BB;
end


x = mids./80;
allM = cell(8,1);
allSEM = cell(8,1);
all95 = cell(8,1);
for con = 1:8
    BBt = BBAll{con};
    mouseN = size(BBt{1,1},1);
    dayN = size(BBt,2);
    
    tempM = cell(L,dayN);
    tempSEM = cell(L,dayN);
    temp95 = cell(L,dayN);
    for d = 1:dayN
        for i = 1:L
            y = BBt{i,d};
            yMean = mean(y,'omitnan');
            ySEM = std(y,'omitnan')/sqrt(mouseN);
            CI95 = tinv([0.025 .975],mouseN-1);
            yCI95 = bsxfun(@times,ySEM,CI95(:));
            
            tempM{i,d} = yMean;
            tempSEM{i,d} = ySEM;
            temp95{i,d} = yCI95;
        end
    end
    allM{con} = tempM;
    allSEM{con} = tempSEM;
    all95{con} = temp95;
end

for b = 1:8
for con = 1:8
    for d= 1:4
    yMean = allM{con}{b,d};
    yCI95 = all95{con}{b,d};
    bm{b}(con,d) = max(max(yMean)) %+yCI95(2,:)));
    end
end
end
for b = 1:8
    maxb(b) = max(max(bm{b}))*1.1;
end


b = 8;
for d = 1:4
figure(d);
for con = [7 8]
yMean = allM{con}{b,d};
yCI95 = all95{con}{b,d};
h = fill([x fliplr(x)],[yMean+yCI95(1,:) fliplr(yMean+yCI95(2,:))],cvcon(con,:));
axis([0 1200 0 maxb(b)])
set(h,'facealpha',.5);
hold on;
end
for con = [7 8]
        yMean = allM{con}{b,d};
plot(x,yMean,'LineWidth',2,'Color',cvcon(con,:))
end
% saveas
%close all
end

% summary by condition
% summary by day for 4 days
w2 = 60*80*20;
B = cell(8,1);
for b = 1:8
    for con = 1:8
        temp = cC1CON{con};
        ttotal = nan(size(temp));
        for d = 1:size(temp,2)
            for m = 1:size(temp,1)
                try
                    tt = temp{m,d};
                    tb = find(tt==b);
                    firstb = tb(1)./(80);
                    totalb = length(tb)./80; % seconds
                    
                    ttotal(m,d) = totalb;
                catch
                end
            end
            B{b}{con} = ttotal;
        end
    end
    
end

con = 7;   
C = cell(8,4); cv = []; mids = []; cv2 = []; w = .6; addPoints = 'true';
for b = 1:8
    mo = (b-1)*5;
    for d = 1:4
        C{b,d} = B{b}{con}(:,d);
        mids = [mids mo+d];
        cv = [cv; 0 0 0];
        cv2 = [cv2; cvcon(con,:)];
    end
end 
        
CD = C'; CD = CD(:); %CD{30}(16) = nan;
oi = myBoxNA(CD,mids,w,cv,addPoints);
mm = 600;%max(combineCells(CD))*1.05;
axis([0 41 0 mm]); hold on;
oim = oi(:,1);
for b = 1:8
    mo = (b-1)*5;
    ds = (1:4)+mo; cid = (1:4)+(b-1)*4;
    scatter(ds,oim(cid),[],cbarbeh(b,:),'filled');
    plot(ds,oim(cid),'Color',cbarbeh(b,:),'LineWidth',3);
end
set(gcf,'Units','Normalized','OuterPosition',[0 0.04 .8 0.5])
%axis off




%% figure 6
% groomingSpecifics.m 
load('data/GROOM_IDX.mat'); % hand-labeled grooms 

%% figure 7
% locomotionSpecifics.m

%% polar scatters - comparisons
jtcolor = [1 0 0; 0 0 1; .8 0 .8; 0 1 1; 1 1 .5; .5 1 0; 0 .8 0];
ai = 2;
% L7-tsc1
figure(1); 
con = 4; allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan');
sz = 500*sz./mz;
polarscatter(allsa(:,1),hx(2:end),sz,'k','filled'); hold on; polarscatter(allsa(:,2),hx(2:end),sz,'k','filled'); hold on;
polarscatter(allsa(:,3),hx(2:end),sz,'k','filled'); hold on; polarscatter(allsa(:,4),hx(2:end),sz,'k','filled'); hold on;
rlim([0 .4]);

con = 5; allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan');
sz = 500*sz./mz;
polarscatter(allsa(:,1),hx(2:end),sz,[.8 .8 .8],'filled'); hold on; polarscatter(allsa(:,2),hx(2:end),sz,[.8 .8 .8],'filled'); hold on;
polarscatter(allsa(:,3),hx(2:end),sz,[.8 .8 .8],'filled'); hold on; polarscatter(allsa(:,4),hx(2:end),sz,[.8 .8 .8],'filled'); hold on;
rlim([0 .4]);

con = 6; allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan'); sz = 500*sz./mz; 
for i = 1:4
polarscatter(allsa(:,i),hx(2:end),sz,jtcolor(i,:),'filled'); hold on; 
end
rlim([0 .4]); 
set(gca,'ThetaTickLabel',[]); set(gca,'RTickLabel',[]);

saveas(gcf,'phaseFigs/L7_tsc1_loc.pdf');
saveas(gcf,'phaseFigs/L7_tsc1_loc.tiff');

% Cntnap2
figure(2);
con = 1; allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan');
sz = 500*sz./mz;
polarscatter(allsa(:,1),hx(2:end),sz,'k','filled'); hold on; polarscatter(allsa(:,2),hx(2:end),sz,'k','filled'); hold on;
polarscatter(allsa(:,3),hx(2:end),sz,'k','filled'); hold on; polarscatter(allsa(:,4),hx(2:end),sz,'k','filled'); hold on;
rlim([0 .4]);

con = 2; allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan');
sz = 500*sz./mz;
polarscatter(allsa(:,1),hx(2:end),sz,[.8 .8 .8],'filled'); hold on; polarscatter(allsa(:,2),hx(2:end),sz,[.8 .8 .8],'filled'); hold on;
polarscatter(allsa(:,3),hx(2:end),sz,[.8 .8 .8],'filled'); hold on; polarscatter(allsa(:,4),hx(2:end),sz,[.8 .8 .8],'filled'); hold on;
rlim([0 .4]);

con = 3; allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan'); sz = 500*sz./mz; 
for i = 1:4
polarscatter(allsa(:,i),hx(2:end),sz,jtcolor(i,:),'filled'); hold on; 
end
rlim([0 .4]);
set(gca,'ThetaTickLabel',[]); set(gca,'RTickLabel',[]);
saveas(gcf,'phaseFigs/Cntnap2_loc.pdf');
saveas(gcf,'phaseFigs/Cntnap2_loc.tiff');

% WT
figure(3);
con = 7; allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan');
sz = 500*sz./mz;
polarscatter(allsa(:,1),hx(2:end),sz,'k','filled'); hold on; polarscatter(allsa(:,2),hx(2:end),sz,'k','filled'); hold on;
polarscatter(allsa(:,3),hx(2:end),sz,'k','filled'); hold on; polarscatter(allsa(:,4),hx(2:end),sz,'k','filled'); hold on;
rlim([0 .4]);

con = 8; allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan');
sz = 500*sz./mz;
polarscatter(allsa(:,1),hx(2:end),sz,[.8 .8 .8],'filled'); hold on; polarscatter(allsa(:,2),hx(2:end),sz,[.8 .8 .8],'filled'); hold on;
polarscatter(allsa(:,3),hx(2:end),sz,[.8 .8 .8],'filled'); hold on; polarscatter(allsa(:,4),hx(2:end),sz,[.8 .8 .8],'filled'); hold on;
rlim([0 .4]);
set(gca,'ThetaTickLabel',[]); set(gca,'RTickLabel',[]);
saveas(gcf,'phaseFigs/WT_loc.pdf');
saveas(gcf,'phaseFigs/WT_loc.tiff');

% con = 3; allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan'); 
% sz = 500*sz./mz; polarscatter(allsa(:,1),hx(2:end),sz,'red','filled'); hold on; polarscatter(allsa(:,2),hx(2:end),sz,'blue','filled'); hold on;
% polarscatter(allsa(:,3),hx(2:end),sz,[152 78 163]./255,'filled'); hold on; polarscatter(allsa(:,4),hx(2:end),sz,'cyan','filled'); hold on;
% rlim([0 .8]);
% 


% comparing alpha bins - angular velocity
jtcolor = [1 0 0; 0 0 1; .8 0 .8; 0 1 1; 1 1 .5; .5 1 0; 0 .8 0];
% L7-tsc1
for con = 1:8
figure(con); 

ai = 1;
allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan');
sz = 500*sz./mz;
polarscatter(allsa(:,1),hx(2:end),sz,'r','filled'); hold on; polarscatter(allsa(:,2),hx(2:end),sz,'r','filled'); hold on;
polarscatter(allsa(:,3),hx(2:end),sz,'r','filled'); hold on; polarscatter(allsa(:,4),hx(2:end),sz,'r','filled'); hold on;
rlim([0 .4]);

ai = 2;
allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan');
sz = 500*sz./mz;
polarscatter(allsa(:,1),hx(2:end),sz,'g','filled'); hold on; polarscatter(allsa(:,2),hx(2:end),sz,'g','filled'); hold on;
polarscatter(allsa(:,3),hx(2:end),sz,'g','filled'); hold on; polarscatter(allsa(:,4),hx(2:end),sz,'g','filled'); hold on;
rlim([0 .4]);

ai = 3;
allsa = infoCON{ai,con}; sz = allsa(:,5); sz(sz==0) = nan; mz = sum(sz,'omitnan');
sz = 500*sz./mz;
polarscatter(allsa(:,1),hx(2:end),sz,'b','filled'); hold on; polarscatter(allsa(:,2),hx(2:end),sz,'b','filled'); hold on;
polarscatter(allsa(:,3),hx(2:end),sz,'b','filled'); hold on; polarscatter(allsa(:,4),hx(2:end),sz,'b','filled'); hold on;
rlim([0 .4]);

%saveas(gcf,['phaseFigs/left_center_right_con_' num2str(con) '.tiff']);
end

%% 






