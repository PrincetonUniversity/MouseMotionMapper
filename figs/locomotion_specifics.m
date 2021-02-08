load('data/allC100_8C_CON.mat','sC1CON'); % 100-> sorting 
load('data/allC100_CON.mat');
load('data/mouseSkeletonSimple_400.mat') % skeleton
load('data/centroids_velocities_CON.mat')

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


CON_Vs = cell(size(CON_V));
for con = 1:8
    ccon = CON_V{con};
    [mN dN] = size(ccon);
    tempv = cell(mN,dN);
    for m = 1:mN
        for d = 1:dN
            vv = ccon{m,d};
            vv2 = vv; vv2(vv>10) = 10;
            tempv{m,d} = vv2;
        end
    end
    CON_Vs{con} = tempv;
end

allL = cell(8,4);
for con = 1:8
    ccon = cC1CON{con};
    vcon = CON_Vs{con};
    [mN dN] = size(ccon);
    onlyL = cell(mN,dN);
    for m = 1:mN
        for d = 1:dN
            cc = ccon{m,d};
            vv = vcon{m,d};
            lv = vv(find(cc==8));
            onlyL{m,d} = lv;
        end
    end
    for d = 1:4
        dol = combineCells(onlyL(:,d));
        allL{con,d} = dol;
    end
end
    
totalTime = zeros(8,4);
for con = 1:8
    for d = 1:4
    totalTime(con,d) = length(combineCells(CON_V{con}(:,d)));
    end
end

lx = 0:.1:.4;
locFracBin = cell(1,8);
for con = 1:8
    temploc = zeros(4,5);
    for d = 1:4
        tempL = allL{con,d}*(80/(1000*1.97));
        temploc(d,:) = hist(tempL,lx)./totalTime(con,d);
    end
    locFracBin{con} = temploc;
end
     
redbarcol = ones(6,3);
redbarcol(:,1) = linspace(1,1,6);
redbarcol(:,2) = linspace(1,0,6);
redbarcol(:,3) = linspace(1,0,6);

lsC2 = [locFracBin{1}; locFracBin{2}; locFracBin{3}];
bx = [1 2 3 4 6 7 8 9 11 12 13 14];
figure(1);
b1 = bar(bx,lsC2,'stacked');
for i = 1:5
    b1(i).FaceColor = 'flat';
    b1(i).CData = redbarcol(i,:);
end
axis([0 15 0 .45]);
%saveas(gcf,['/Users/ugne/Dropbox/mouseMaster/loc_C2.pdf'])
axis off
%saveas(gcf,['/Users/ugne/Dropbox/mouseMASTER/loc_C2_NA.pdf'])
%close(gcf)

lsL7 = [locFracBin{4}; locFracBin{5}; locFracBin{6}];
bx = [1 2 3 4 6 7 8 9 11 12 13 14];
figure(2);
b1 = bar(bx,lsL7,'stacked');
for i = 1:5
    b1(i).FaceColor = 'flat';
    b1(i).CData = redbarcol(i,:);
end
axis([0 15 0 .45]);
%saveas(gcf,['/Users/ugne/Dropbox/mouseMaster/loc_L7.pdf'])
axis off
%saveas(gcf,['/Users/ugne/Dropbox/mouseMASTER/loc_L7_NA.pdf'])
%close(gcf)


lsWT = [locFracBin{7}; locFracBin{8}];
bx = [1 2 3 4 6 7 8 9];
figure(3);
b1 = bar(bx,lsWT,'stacked');
for i = 1:5
    b1(i).FaceColor = 'flat';
    b1(i).CData = redbarcol(i,:);
end
axis([0 10 0 .45]);
%saveas(gcf,['/Users/ugne/Dropbox/mouseMaster/loc_WT.pdf'])
axis off
%saveas(gcf,['/Users/ugne/Dropbox/mouseMASTER/loc_WT_NA.pdf'])
%close(gcf)

% locomotion phase plots

% joints are saved in variable jointsCON ->
% this is a 1x8 cell array corresponding to each condition (cntnap2(1,2,3),
% L7-tsc1(4,5,6), Black6 males (7) and Black6 females (8)
% each cell holds a [numMouse x numDay] cell array corresponding to mouse
% and day for each trial in that condition
% jointsCON{condition}{mouseNum,dayNum} -> 18x2xt array of joint positions
% corresponding to skeleton

% STEP 1: 
juse = [7 8 14 15 16 17 18];
jusexy = [juse juse+18];
% just locomotion sorting
MCLCON = cell(1,8);
for con = 1:8
[mN dN] = size(jointsCON{con});
MCL = cell(mN, dN);
for m = 1:mN
    m
    for d = 1:dN
        try
        ja20 = jointsCON{con}{m,d};
        r100 = cC1CON{con}{m,d};
        %r100 = sortid8(c100);
        % sort runs through the center and edges
        cc = bwconncomp(r100==8); n = length(r100);
        
        cl = zeros(length(cc.PixelIdxList),1);
        for i = 1:length(cl)
            cl(i) = length(cc.PixelIdxList{i});
        end
        
        c2.PixelIdxList = cc.PixelIdxList(cl>=80);
        lc2 = length(c2.PixelIdxList);
        p1Center = findP1Center(ja20);
        
        i = 9; jt1 = squeeze(ja20(i,:,:)); jt1(1,:) = smooth(medfilt1(jt1(1,:),4),40); jt1(2,:) = smooth(medfilt1(jt1(2,:),4),40);
        i = 16; jt2 = squeeze(ja20(i,:,:)); jt2(1,:) = smooth(medfilt1(jt2(1,:),4),40); jt2(2,:) = smooth(medfilt1(jt2(2,:),4),40);
        centpt = zeros(n,2); heada = zeros(n,1);
        centpt(:,1) = (jt1(1,:)+jt2(1,:))./2; centpt(:,2) = (jt1(2,:)+jt2(2,:))./2;
        dir = jt2-jt1;
        for i = 1:n
            heada(i) = atan2d(dir(2,i),dir(1,i));
        end
        h2a = smooth(unwrap(heada),20);
        headav = [0; diff(h2a)]; 
        
        vt = [0 0; diff(centpt)]; vt = sqrt(sum(vt.^2,2)); cent = vt;
        vt(vt>10) = nan; vt = fillmissing(vt,'linear'); vt = smooth(vt,20);
        linfo = cell(lc2,7);
        clear pt ptc ptmid dtrav dd ts cvel avel acfrac
        
        for i = 1:lc2
            pt = ja20(:,:,c2.PixelIdxList{i});
            ptc = p1Center(jusexy,c2.PixelIdxList{i});
            cvel = vt(c2.PixelIdxList{i});
            linfo{i,1} = ptc;
            linfo{i,2} = reshape(pt,[36 size(pt,3)]);
            linfo{i,3} = cvel;
            linfo{i,4} = centpt(c2.PixelIdxList{i},:);
            linfo{i,5} = heada(c2.PixelIdxList{i});
            linfo{i,6} = headav(c2.PixelIdxList{i});
            linfo{i,7} = c2.PixelIdxList{i};
        end
        MCL{m,d} = linfo;
        catch
        end
    end
end
MCLCON{con} = MCL;
end
% save - all packaged info per snippet
% save('MCLCON_allinfo.mat','MCLCON','-v7.3');
load('data/MCLCON_allinfo.mat','MCLCON');

% STEP2
% have locoCONall - now collect data for each stride in each snippet

% cv - centroid velocity
% jv - joint velocities
% trace - all joint traces - real space
% rose - phase angle of each paw at timepoint captured
% idx - frame index for timepoints used
% day - recording day number
% mouse - mouse number within condition
% place - centroid coord @ phase 0
% placeA - centroid coordinates for entire snippet
% mtj - centered/rotated coordinates of joints used (juse)
% heada - heading
% shape - path traveled
% turn - heading angular velocity

locoCONall = cell(1,8);
for con = 1:8
    locoCON.cv = [];
    locoCON.jv = [];
    locoCON.trace = [];
    locoCON.rose = [];
    locoCON.idx = [];
    locoCON.day = [];
    locoCON.mouse = [];
    locoCON.place = [];
    locoCON.placeA = [];
    locoCON.mtj = [];
    locoCON.heada = [];
    locoCON.shape = [];
    locoCON.turn = [];
    ccount = 1;
    
    mclcon = MCLCON{con};
    [mN dN] = size(mclcon);
    for m = 1:mN
        fprintf(1,['Processing mouse # ' num2str(m) '\n']);
        for d = 1:dN
            mt = mclcon{m,d};
            jN = size(mt,1);
            for j = 1:jN
                
                mtj = mt{j,1};
                
                mtji = zeros(size(mtj,1),size(mtj,2)*100);
                for x = 1:size(mtj,1)
                mtji(x,:) = interp(mtj(x,:),100); 
                end
                
                mt5 = mt{j,2}; vt5 = zeros(7,size(mt5,2));
                vt5(1,:) = [0 sqrt(diff(squeeze(mt5(7,:))).^2 + diff(squeeze(mt5(7+18,:))).^2)];
                vt5(2,:) = [0 sqrt(diff(squeeze(mt5(8,:))).^2 + diff(squeeze(mt5(8+18,:))).^2)];
                vt5(3,:) = [0 sqrt(diff(squeeze(mt5(14,:))).^2 + diff(squeeze(mt5(14+18,:))).^2)];
                vt5(4,:) = [0 sqrt(diff(squeeze(mt5(15,:))).^2 + diff(squeeze(mt5(15+18,:))).^2)];
                vt5(5,:) = [0 sqrt(diff(squeeze(mt5(16,:))).^2 + diff(squeeze(mt5(16+18,:))).^2)];
                vt5(6,:) = [0 sqrt(diff(squeeze(mt5(17,:))).^2 + diff(squeeze(mt5(17+18,:))).^2)];
                vt5(7,:) = [0 sqrt(diff(squeeze(mt5(18,:))).^2 + diff(squeeze(mt5(18+18,:))).^2)];
                
                ztj = zscore(mtji')';
                [p1 id1] = findpeaks(ztj(8,:),'MinPeakProminence',.2,'MinPeakDistance',900);
                [p2 id2] = findpeaks(ztj(9,:),'MinPeakProminence',.2,'MinPeakDistance',900);
                [p3 id3] = findpeaks(ztj(10,:),'MinPeakProminence',.2,'MinPeakDistance',900);
                [p4 id4] = findpeaks(ztj(11,:),'MinPeakProminence',.2,'MinPeakDistance',900);
                
                zn = size(ztj,2);
                phasetemp = nan(4,zn);

                
                % rewrite phase-getting to check between 0-2pi for all paws
                sl = size(mtji,2);
                nid1 = length(id1)-1;
                for i = 1:nid1
                    lt = length(id1(i):id1(i+1));
                    tt = linspace(0,360,lt);
                    phasetemp(1,id1(i):id1(i+1)) = tt;
                end
                
                nid2 = length(id2)-1;
                for i = 1:nid2
                    lt = length(id2(i):id2(i+1));
                    tt = linspace(0,360,lt);
                    phasetemp(2,id2(i):id2(i+1)) = tt;
                end
                
                nid3 = length(id3)-1;
                for i = 1:nid3
                    lt = length(id3(i):id3(i+1));
                    tt = linspace(0,360,lt);
                    phasetemp(3,id3(i):id3(i+1)) = tt;
                end
                
                nid4 = length(id4)-1;
                for i = 1:nid4
                    lt = length(id4(i):id4(i+1));
                    tt = linspace(0,360,lt);
                    phasetemp(4,id4(i):id4(i+1)) = tt;
                end
                
                p11 = find(phasetemp(1,:)==0);
                p12 = find(phasetemp(1,:)==360);
                use1 = [p11 p12];
                
 
                % pn1 = round(id1(find(use1==1))./100);
                pcon = phasetemp(:,use1);
                pn1 = round(use1./100);
                locoCON.day{ccount} = d;
                locoCON.mouse{ccount} = m;
                locoCON.cv{ccount} = mt{j,3};
                locoCON.jv{ccount} = vt5;
                locoCON.trace{ccount} = mt5;
                locoCON.rose{ccount} = pcon;
                locoCON.idx{ccount} = pn1;
                locoCON.place{ccount} = mt{j,4}(pn1,:);
                locoCON.placeA{ccount} = mt{j,4};
                locoCON.mtj{ccount} = mtj;
                locoCON.heada{ccount} = mt{j,5};
                locoCON.turn{ccount} = mt{j,6};
                locoCON.shape{ccount} = mt{j,4}-mt{j,4}(1,:);
                ccount = ccount+1;
            end
        end
    end
    locoCONall{con} = locoCON;
end

%save('julyFigs/locoCONall.mat','locoCONall','-v7.3');
% STEP 3
% sort through snippets to find info and phase matches for matched time
% points

stepsCON = cell(1,8); phaseCON = cell(1,8);
for con = 1:8
    tcon = locoCONall{con};
    ni = length(tcon.cv);
    paw1 = [];
    paw2 = [];
    paw3 = [];
    paw4 = [];
    day = [];
    mouse = [];
    vv = [];
    place = [];
    headv = [];
    ssnum = [];
    % check here for offsets
    for i = 1:ni
        cv = tcon.cv{i}.*(80/(1.97*1000)); jv = tcon.jv{i}; trace = tcon.trace{i}; rose = tcon.rose{i};
        idx = tcon.idx{i}+1; d = tcon.day{i}; cent = tcon.place{i}; m = tcon.mouse{i}; hv = tcon.turn{i}(idx);
        
        ss = size(rose,2);
        paw1 = [paw1; deg2rad(rose(1,:))'];
        paw2 = [paw2; deg2rad(rose(2,:))'];
        paw3 = [paw3; deg2rad(rose(3,:))'];
        paw4 = [paw4; deg2rad(rose(4,:))'];
        day = [day; repmat(d,[ss 1])];
        mouse = [mouse; repmat(m,[ss 1])];
        vv = [vv; cv(idx)];
        place = [place; cent];
        headv = [headv; hv];
        ssnum = [ssnum; [1:ss]'];
    end
    infoCON = [mouse day ssnum vv headv place];
    specs =  [paw1 paw2 paw3 paw4];
    stepsCON{con} = infoCON; phaseCON{con} = specs;
end


% specify bin sizes here for centroid velocity and angular velocity
cnum = 8; % number of conditions
abins = [-3 -1 1 3]; al = length(abins)-1;
hx = linspace(0,.4,21);
infoCON = cell(al,cnum);
roseCON = cell(al,cnum);

for con = 1:cnum
    consteps = stepsCON{con}; 
    conphase = phaseCON{con};

        for ai = 1:al
            tempa = consteps(:,5);
            ax = find(tempa>abins(ai) & tempa<abins(ai+1));
            vsort = consteps(ax,4);
            
            tsa = nan(20,4); nsa = zeros(20,1); vsa = nan(20,4);
            for j = 1:20
                useic = find(hx(j)<=vsort & vsort<=hx(j+1));
                temps = consteps(ax(useic),:);
                tphase = conphase(ax(useic),:);
                
                tphaseu = sum(~isnan(tphase),2);
                utp = find(tphaseu==4);
                
                tphase1 = tphase(utp,:);
                tsa(j,:) = median(tphase1,1); nsa(j) = size(tphase1,1); vsa(j,:) = std(tphase1,1);
            end
            
            allsa = [tsa nsa vsa];
            infoCON{ai,con} = allsa;
            roseCON{ai,con} = tsa;
        end
        
end

haicon = zeros(8,61);
for con = 1:8
    consteps = stepsCON{con}; 
    conphase = phaseCON{con};
    aii = consteps(:,5);
    haicon(con,:) = hist(aii,-3:.1:3)./length(aii);
end

for i = 1:8
    plot(haicon(i,:),'Color',cvcon(i,:),'LineWidth',2); hold on;
end


% plot all angle bins for black6 males
con = 7;
%infoCON = stepsCON{con};
for ai = 1:3
subplot(3,1,ai);
allsa = infoCON{ai,con};
sz = allsa(:,5); sz(sz==0) = nan;
mz = sum(sz,'omitnan');

sz = 300*sz./mz;
polarscatter(allsa(:,1),hx(2:end),sz,'red','filled'); hold on;
polarscatter(allsa(:,2),hx(2:end),sz,'blue','filled'); hold on;
polarscatter(allsa(:,3),hx(2:end),sz,[152 78 163]./255,'filled'); hold on;
polarscatter(allsa(:,4),hx(2:end),sz,'cyan','filled'); hold on;
rlim([0 .4])
end


% only use middle bin for comparing across groups for now
ai = 2;
