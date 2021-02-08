%% wild type grooming patterns
%load('manuscriptMOUSE/jointsCON.mat'); % joint positions in real space
load('data/GROOM_IDX.mat'); % hand-labeled grooms 
load('data/allC100_8C_CON.mat','sC1CON'); % 100-> sorting 
load('data/mouseSkeletonSimple_400.mat') % skeleton


G = cell(size(jointsCON));

for con = 1:8
    [mN dN] = size(jointsCON{con});
    MCL = cell(mN, dN);

    for d = 1:dN
        for m = 1:mN
            m
            try
                ja20 = jointsCON{con}{m,d}; n = length(ja20);
                c100 = cC1CON{con}{m,d};
                f100 = fC1CON{con}{m,d};
                
                cc = bwconncomp(c100==2);
                cl = zeros(length(cc.PixelIdxList),1);
                for i = 1:length(cl)
                    cl(i) = length(cc.PixelIdxList{i});
                end
                c2.PixelIdxList = cc.PixelIdxList(cl>=20);
                lc2 = length(c2.PixelIdxList);
                
                tv = zeros(18,n);
                for i = 1:18
                    tj = squeeze(ja20(i,:,:))';
                    vv = [0 0; diff(tj)];
                    tv(i,:) = sqrt(sum(vv.^2,2));
                end
                tvn = tv.*(80/(1000*1.97));
                tvn(tvn>1) = 1;
                
                ccid = zeros(size(c2.PixelIdxList,2),3);
                for ji = 1:size(c2.PixelIdxList,2)
                    ccid(ji,1) = c2.PixelIdxList{ji}(1);
                    ccid(ji,2) = c2.PixelIdxList{ji}(end);
                    ccid(ji,3) = c2.PixelIdxList{ji}(end)-c2.PixelIdxList{ji}(1);
                end
                
                gg.pos = [];
                gg.len = [];
                gg.p3 = [];
                gg.usage = [];
                gg.ccid = [];
                gg.x = [];
                for ii = 1:size(ccid,1)
                    %figure(ii);
                    tvns = sum(tvn(:,ccid(ii,1):ccid(ii,2)),2)./(ccid(ii,2)-ccid(ii,1));
                    %mean posture over groom snippet,
                    %bubbles over limbs - usage
                    posS = ja20(:,:,ccid(ii,1):ccid(ii,2));
                    posm = zeros(18,2);
                    for i = 1:18
                        ps = squeeze(posS(i,:,:))';
                        posm(i,:) = median(ps);
                    end
                    tailP = posm(16,:); centerP = posm(9,:);
                    posm2 = posm-repmat(tailP,[18 1]);
                    [tx rx] = cart2pol(posm2(:,1),posm2(:,2));
                    xtemp = tailP-centerP;
                    rotdir1 = deg2rad(rem(atan2d(xtemp(1),xtemp(2))+360,360));
                    
                    [pos3(:,1),pos3(:,2)] = pol2cart(tx+rotdir1,rx);
                    
                    for i = 1:length(segments.joints)
                        j1 = segments.joints_idx{i}(1);
                        j2 = segments.joints_idx{i}(2);
                        p = [pos3(j1,:); pos3(j2,:)];
                        %plot(p(:,1),p(:,2)); hold on;
                    end
%                     for i = 1:18
%                         scatter(pos3(i,1),pos3(i,2),tvns(i)*100,'r','filled'); hold on;
%                     end
%                     for i = 1:18
%                         scatter(pos3(i,1),pos3(i,2),'k'); hold on
%                         %text(posm(i,1)+2,posm(i,2)+2,joints.name{i});
%                     end
%                     axis equal ij
%                     
                    gg.pos{ii} = posm;
                    gg.x{ii} = tailP;
                    gg.len{ii} = ccid(ii,3);
                    gg.p3{ii} = pos3;
                    gg.usage{ii} = tvns;
                    gg.ccid{ii} = ccid(ii,:);
                    
                end
                
                G{con}{m,d} = gg;
                
            catch
            end
        end
    end
end




% total entries into grooming
% lengths of bouts
GPX = cell(1,8);
GPT = cell(1,8);
GPG = cell(1,8);
for con = 1:8
    con
    [mN dN] = size(G{con});
    for m = 1:mN
        for d = 1:dN
            
            try
                gg = G{con}{m,d};
                lx = length(gg.x); gpx = zeros(lx,2); gpt = zeros(lx,2);
                for i = 1:lx
                    gpx(i,:) = gg.x{i};
                    gpt(i,:) = gg.ccid{i}([1 3]);
                end
                tlink = zeros(lx-1,1); dlink = zeros(lx-1,1);
                for ji = 1:lx-1
                    tlink(ji) = gg.ccid{ji+1}(1)-gg.ccid{ji}(2);
                    dlink(ji) = pdist2(gpx(ji,:),gpx(ji+1,:));
                end
                groupings = zeros(lx,1);
                groupings(1) = 1; gcurr = 1;
                
                for ji = 2:lx
                    if tlink(ji-1) < 800 & dlink(ji-1) < 10
                        groupings(ji) = gcurr;
                    else
                        gcurr = gcurr+1;
                        groupings(ji) = gcurr;
                    end
                end
                GPX{con}{m,d} = gpx;
                GPT{con}{m,d} = gpt;
                GPG{con}{m,d} = groupings;
            catch
            end
        end
    end
end

wall1(:,1) = 120:920; wall1(:,2) = 250*ones(size(120:920));
wall2(:,1) = 120*ones(size(250:1050)); wall2(:,2) = 250:1050;
wall3(:,1) = 120:920; wall3(:,2) = 1050*ones(size(120:920));
wall4(:,1) = 920*ones(size(250:1050)); wall4(:,2) = 250:1050;


for con = 1:8
    figure(con)
    for m = 1:10
        try
        for d = 1:4
            subplot(10,4,(m-1)*4+d);
            gpx = GPX{con}{m,d}; gpt = GPT{con}{m,d};
        groupings = GPG{con}{m,d};
        scatter(gpx(:,1),gpx(:,2),gpt(:,2),gpt(:,1)); hold on;
        colormap(cmap3); axis equal off
        plot(wall1(:,1),wall1(:,2),'k');
        plot(wall2(:,1),wall2(:,2),'k');
        plot(wall3(:,1),wall3(:,2),'k');
        plot(wall4(:,1),wall4(:,2),'k');
        end
        catch
        end
    end
    set(gcf,'Units','Normalized','OuterPosition',[0 0.04 .5 0.96])
    %saveas(gcf,['augustFigs/groomEx/con_' num2str(con) '.tiff']);

end


con = 1;
for m = 1:14
    figure(m);
    for d = 1:4
        subplot(1,4,d);
        gpx = GPX{con}{m,d}; gpt = GPT{con}{m,d};
        groupings = GPG{con}{m,d};
        scatter(gpx(:,1),gpx(:,2),gpt(:,2),gpt(:,1)); hold on;
        colormap(cmap3); axis equal off
        plot(wall1(:,1),wall1(:,2),'k');
        plot(wall2(:,1),wall2(:,2),'k');
        plot(wall3(:,1),wall3(:,2),'k');
        plot(wall4(:,1),wall4(:,2),'k');
    end
end

con = 2;
for m = 1:15
    figure(m);
    for d = 1:4
        subplot(1,4,d);
        gpx = GPX{con}{m,d}; gpt = GPT{con}{m,d};
        groupings = GPG{con}{m,d};
        scatter(gpx(:,1),gpx(:,2),gpt(:,2),gpt(:,1)); hold on;
        colormap(cmap3); axis equal off
        plot(wall1(:,1),wall1(:,2),'k');
        plot(wall2(:,1),wall2(:,2),'k');
        plot(wall3(:,1),wall3(:,2),'k');
        plot(wall4(:,1),wall4(:,2),'k');
    end
end

con = 3;
for m = 1:10
    figure(m);
    for d = 1:4
        subplot(1,4,d);
        gpx = GPX{con}{m,d}; gpt = GPT{con}{m,d};
        groupings = GPG{con}{m,d};
        scatter(gpx(:,1),gpx(:,2),gpt(:,2),gpt(:,1)); hold on;
        colormap(cmap3); axis equal off
        plot(wall1(:,1),wall1(:,2),'k');
        plot(wall2(:,1),wall2(:,2),'k');
        plot(wall3(:,1),wall3(:,2),'k');
        plot(wall4(:,1),wall4(:,2),'k');
    end
end


con = 4;
for m = 1:17
    figure(m);
    for d = 1:4
        subplot(1,4,d);
        gpx = GPX{con}{m,d}; gpt = GPT{con}{m,d};
        groupings = GPG{con}{m,d};
        scatter(gpx(:,1),gpx(:,2),gpt(:,2),gpt(:,1)); hold on;
        colormap(cmap3); axis equal off
        plot(wall1(:,1),wall1(:,2),'k');
        plot(wall2(:,1),wall2(:,2),'k');
        plot(wall3(:,1),wall3(:,2),'k');
        plot(wall4(:,1),wall4(:,2),'k');
    end
end

con = 5;
for m = 1:17
    figure(m);
    for d = 1:4
        subplot(1,4,d);
        gpx = GPX{con}{m,d}; gpt = GPT{con}{m,d};
        groupings = GPG{con}{m,d};
        scatter(gpx(:,1),gpx(:,2),gpt(:,2),gpt(:,1)); hold on;
        colormap(cmap3); axis equal off
        plot(wall1(:,1),wall1(:,2),'k');
        plot(wall2(:,1),wall2(:,2),'k');
        plot(wall3(:,1),wall3(:,2),'k');
        plot(wall4(:,1),wall4(:,2),'k');
    end
end

con = 6;
for m = 1:9
    figure(m);
    for d = 1:4
        subplot(1,4,d);
        gpx = GPX{con}{m,d}; gpt = GPT{con}{m,d};
        groupings = GPG{con}{m,d};
        scatter(gpx(:,1),gpx(:,2),gpt(:,2),gpt(:,1)); hold on;
        colormap(cmap3); axis equal off
        plot(wall1(:,1),wall1(:,2),'k');
        plot(wall2(:,1),wall2(:,2),'k');
        plot(wall3(:,1),wall3(:,2),'k');
        plot(wall4(:,1),wall4(:,2),'k');
    end
end


con = 7;
for m = 1:60
    figure(m);
    for d = 1:4
        subplot(1,4,d);
        gpx = GPX{con}{m,d}; gpt = GPT{con}{m,d};
        groupings = GPG{con}{m,d};
        scatter(gpx(:,1),gpx(:,2),gpt(:,2),gpt(:,1)); hold on;
        colormap(cmap3); axis equal off
        plot(wall1(:,1),wall1(:,2),'k');
        plot(wall2(:,1),wall2(:,2),'k');
        plot(wall3(:,1),wall3(:,2),'k');
        plot(wall4(:,1),wall4(:,2),'k');
    end
end

con = 8;
for m = 1:20
    figure(m);
    for d = 1:4
        try
        subplot(1,4,d);
        gpx = GPX{con}{m,d}; gpt = GPT{con}{m,d};
        groupings = GPG{con}{m,d};
        scatter(gpx(:,1),gpx(:,2),gpt(:,2),gpt(:,1)); hold on;
        colormap(cmap3); axis equal off
        plot(wall1(:,1),wall1(:,2),'k');
        plot(wall2(:,1),wall2(:,2),'k');
        plot(wall3(:,1),wall3(:,2),'k');
        plot(wall4(:,1),wall4(:,2),'k');
        catch
        end
    end
end


% number groups and length of bouts, fraction of actual groom
% gstat: startTime, totalTime, groomTime, groomFrac
gsummary = cell(1,8);
gsummary2 = cell(1,8);
for con = 1:8
    [mN dN] = size(GPX{con});
    for m = 1:mN
        for d = 1:dN
            gpx = GPX{con}{m,d};
            gpt = GPT{con}{m,d};
            groupings = GPG{con}{m,d};
            ng = max(groupings);
            
            gstat = zeros(ng,4);
            for n = 1:ng
                fx = find(groupings == n);
                tt = gpt(fx(end),1) - gpt(fx(1),1) + gpt(fx(end),2);
                grooml = sum(gpt(fx,2));
                gfrac = grooml/tt; startt = gpt(fx(1),1);
                gstat(n,:) = [startt tt grooml gfrac];
            end
            gstat2 = gstat((gstat(:,2)>80),:);
            gsummary{con}{m,d} = gstat;
            gsummary2{con}{m,d} = gstat2;
        end
    end
end


% three-colored bars over time
% 0. other 
% 1. between
% 2. groom 
C = cell(1,8);
for con = 1:8
    ccon = sC1CON{con};
    [mN dN] = size(ccon);
    n = 95998;
    temp = zeros(mN,dN,n);
    for d= 1:dN
        for m = 1:mN
            gs = gsummary{con}{m,d};
            ctemp = ccon{m,d};
            ge = gs(:,1)+gs(:,2);
            ste = [gs(:,1) ge];
            
            for i = 1:size(ste,1)
                temp(m,d,ste(i,1):ste(i,2)) = 1;
            end
            temp(m,d,find(ctemp==2)) = 2;
        end
    end
    C{con} = temp;
end

for con = 1:8
    figure(con)
for d = 1:4
subplot(4,1,d)
imagesc(squeeze(C{con}(:,d,:)))
colormap(cmap1)
end
end

offd = [0 100000 200000 300000];
for con = 1:8
    figure(con);
    [mN dN] = size(gsummary2{con});
    
for d = 1:4
    od = offd(d);
for m = 1:mN
gs = gsummary2{con}{m,d};
scatter(gs(:,1)+od,gs(:,2),[],'k','.'); hold on;
%scatter(gs(:,1)+od,gs(:,3),[],'r','filled');
end
end
axis([0 400000 0 9600])
end


dcol = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
for con = 1:8
    figure(con);
    [mN dN] = size(gsummary2{con});
    
for d = 1:4
    od = offd(d);
for m = 1:mN
gs = gsummary2{con}{m,d};
scatter(gs(:,1)+od,gs(:,2),[],dcol(d,:),'.'); hold on;
%scatter(gs(:,1)+od,gs(:,3),[],'r','filled');
end
end
axis([0 400000 0 9600])
end



% median groom length at time bins

for con = 1:8
    %figure(con);
    [mN dN] = size(gsummary2{con});
    
for d = 1:dN
gsc = [];
for m = 1:mN
gs = gsummary2{con}{m,d};
gsc = [gsc; gs];
%scatter(gs(:,1)+od,gs(:,2),[],dcol(d,:),'.'); hold on;
%scatter(gs(:,1)+od,gs(:,3),[],'r','filled');
end
gscd{d} = gsc;
end

cgscd{con} = gscd;
end



tj1 = squeeze(tempj(5,:,:))';


% entries per individual 
% length per entry
ls = zeros(8,4);
lt = zeros(8,4);
ll = zeros(8,4);
lf = zeros(8,4);
for con = 1:8
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    
    for d = 1:dN
        ls(con,d) = size(temp{1,d},1)./mN;
        lt(con,d) = median(temp{1,d}(:,1))./(80*60); %minutes
        ll(con,d) = median(temp{1,d}(:,2))./(80); %seconds
        lf(con,d) = median(temp{1,d}(:,4)); 
    end
end



% # of entries, time spent

cvcon = [188 128 189; 103 169 207; 2 129 138; ...
    253 141 60; 227 26 28; 128 0 38; 0 0 0; 248 129 191]./256;

figure(1);
subplot(1,3,1);
for con = 1:3
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,ls(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,ls(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(ls))]);
subplot(1,3,2);
for con = 4:6
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,ls(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,ls(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(ls))]);
subplot(1,3,3);
for con = 7:8
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,ls(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,ls(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(ls))]);


figure(2);
subplot(1,3,1);
for con = 1:3
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,lt(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,lt(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(lt))]);
subplot(1,3,2);
for con = 4:6
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,lt(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,lt(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(lt))]);
subplot(1,3,3);
for con = 7:8
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,lt(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,lt(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(lt))]);




figure(3);
subplot(1,3,1);
for con = 1:3
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,ll(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,ll(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(ll))]);
subplot(1,3,2);
for con = 4:6
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,ll(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,ll(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(ll))]);
subplot(1,3,3);
for con = 7:8
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,ll(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,ll(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(ll))]);




figure(4);
subplot(1,3,1);
for con = 1:3
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,lf(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,lf(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(lf))]);
subplot(1,3,2);
for con = 4:6
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,lf(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,lf(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(lf))]);
subplot(1,3,3);
for con = 7:8
    temp = cgscd{con};
    [mN dN] = size(sC1CON{con});
    scatter(1:dN,lf(con,1:dN),[],cvcon(con,:),'filled'); hold on;
    plot(1:dN,lf(con,1:dN),'LineWidth',.5,'Color',cvcon(con,:));
end
axis([0 6 0 1.1*max(max(lf))]);


conZone = zeros(8,3);



for con = 1:8
    cd = CON_P{con};
    tempc = CON_P{con};
    tempb = cC1CON{con};
    alle = combineCells(tempc(:)');
    allb = combineCells(tempb(:)');
    bgroom = find(allb==2);
    cgroom = alle(bgroom,:);
    n = histcounts2(cgroom(:,1),cgroom(:,2),0:1200,0:1200);
    nn = n./(sum(sum(n)));
    z1 = nn.*(m1);
    z2 = nn.*(m2);
    z3 = nn.*(m3);
    conZone(con,1) = sum(sum(z1));
    conZone(con,2) = sum(sum(z2));
    conZone(con,3) = sum(sum(z3));
end

for i = 1:8
    bar(i,conZone(i,3),'FaceColor',cvcon(i,:)); hold on;
end
axis off



day1C = cell(1,8);


%%
juse = 1:18;
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
        ja20 = jointsCON{con}{m,d}; n = length(ja20);
        c100 = cC1CON{con}{m,d};
        f100 = fC1CON{con}{m,d};
        
        cc = bwconncomp(c100==2);
        cl = zeros(length(cc.PixelIdxList),1);
        for i = 1:length(cl)
            cl(i) = length(cc.PixelIdxList{i});
        end
        %c2.PixelIdxList = cc.PixelIdxList(cl>=20);
        %lc2 = length(c2.PixelIdxList);
        % p1Center = findP1Center(ja20);
        
        
        
        
        ccid = zeros(size(cc.PixelIdxList,2),2);
        for ji = 1:size(cc.PixelIdxList,2)
            ccid(ji,1) = cc.PixelIdxList{ji}(1);
            ccid(ji,2) = cc.PixelIdxList{ji}(end);
        end
        
        tlink = zeros(size(cc.PixelIdxList,2)-1,1); dlink = zeros(size(cc.PixelIdxList,2)-1,1);
        for ji = 1:size(cc.PixelIdxList,2)-1
           tlink(ji) = ccid(ji+1,1)-ccid(ji,2); 
           dlink(ji) = pdist2(centpt(cc.PixelIdxList{ji+1}(1),:),centpt(cc.PixelIdxList{ji}(end),:));
        end
        groupings = zeros(length(cc.PixelIdxList),1);
        groupings(1) = 1; gcurr = 1;
        for ji = 2:length(cc.PixelIdxList)
            if tlink(ji-1) < 400 & dlink(ji-1) < 10
                groupings(ji) = gcurr;
            else
                gcurr = gcurr+1;
                groupings(ji) = gcurr;
            end
        end
        
        ngroups = gcurr; groomid = zeros(length(ngroups),2);
        for ji = 1:ngroups
            jx = find(groupings==ji);
            groomid(ji,1) = cc.PixelIdxList{jx(1)}(1);
            groomid(ji,2) = cc.PixelIdxList{jx(end)}(end);
            groomid(ji,3) = groomid(ji,2)-groomid(ji,1);
        end
        
        groomid2 = groomid(find(groomid(:,3)>80),:);
        
        
        lc2 = length(groomid2);
        linfo = cell(lc2,8);
        clear pt ptc ptmid dtrav dd ts cvel avel acfrac
        for i = 1:lc2
            currid = groomid2(i,1):groomid2(i,2);
            pt = ja20(:,:,currid);
            %ptc = p1Center(jusexy,currid);
            cvel = vt(currid);
            %linfo{i,1} = ptc;
            linfo{i,1} = pt;
            linfo{i,2} = cvel;
            linfo{i,3} = centpt(currid,:);
            linfo{i,4} = heada(currid);
            linfo{i,5} = headav(currid);
            linfo{i,6} = [currid(1),currid(end)];
            linfo{i,7} = fC1CON{con}{m,d}(currid);
        end
        MCL{m,d} = linfo;
        catch
        end
    end
end
MCLCON{con} = MCL;
end


%%

con = 7;
MCL = MCLCON{con};

gcond = cell(8,4);

for con = 1:8
    MCL = MCLCON{con};
    for d = 1:4
        temp = MCL(:,d);
        tt = [];
        for m = 1:length(temp)
            try
            t1 = temp{m}(:,7);
            tt = [tt; t1];
            catch
            end
        end
        gcond{con,d} = tt;
    end
end

groom3 = cell(8,4); groom4 = cell(8,4);
for con = 1:8
    for d = 1:4
        gt = gcond{con,d};
        ni = length(gt);
        nn = zeros(1,ni);
        for i = 1:ni
            nn(i) = length(gt{i});
        end
        [ss sn] = sort(nn);
        sortedg = gt(sn);
        xa = max(ss);
        newm = nan(ni,xa);
        for i = 1:ni
            newm(i,1:length(sortedg{i})) = sortedg{i};
        end
        newm2 = nan(ni,xa);
        for i = 1:ni
            newm2(i,end-length(sortedg{i})+1:end) = sortedg{i};
        end
        newm3 = nan(size(newm2));
        newm3(~isnan(newm2))=0;
        for i = 1:100
            newm3(find(newm2==i)) = groomids(i);
        end
        groom3{con,d} = newm3;
        newm4 = nan(size(newm3,1),2400);
        try
            newm4(:,end-size(newm3,2)+1:end) = newm3;
        catch
            newm4 = newm3(:,end-2399:end);
        end
        groom4{con,d} = newm4;
    end
end

% nall = newm2(:); nnall = nall(~isnan(nall));



load('/Users/ugne/Dropbox/manuscriptMOUSE/info_8con_8b.mat')
groomid = find(sortid8==2);

notClean = [7 17 26 99];
body = [3 33 36 42 44 50 56 58 61 84 92 96];
nosePaws = [10 15 38 98];
quickSmall = [40 65];
footScratch = [87];




groomids = zeros(1,100);
groomids(body) = 1;
groomids(nosePaws) = 2;
groomids(quickSmall) = 3;
groomids(footScratch) = 4;
groomids(notClean) = 5;
% notClean = 5;


d = 2;
for con = 1:8
    subplot(8,1,con);
    imagesc(groom4{con,d}); colormap(cmapLL); caxis([-.1 5]);
end

gcon = cell(1,8);
for con = 1:8
    [mN dN] = size(fC1CON{con});
    gg = zeros(mN,5,dN);
    for m = 1:mN
        for d = 1:dN
            temp = fC1CON{con}{m,d};
            groomtemp = zeros(size(temp));
            for i = 1:100
                groomtemp(temp==i) = groomids(i);
            end
            gidx = zeros(1,5);
            for i = 1:5
                gidx(i) = sum(groomtemp==i)./length(groomtemp);
            end
            gg(m,:,d) = gidx;
        end
    end
    gcon{con} = gg;
end

%% JOSH HERE 
save('groomMatrix.mat','gcon')


d = 1;
for con = 1:8
    subplot(3,3,con); imagesc(squeeze(gcon{con}(:,:,d)));
    caxis([0 .15])
end

for d = 1:4
    figure(d)
    for con = 1:8
        subplot(8,1,con)
        bar(gcon{con}(:,1:5,d),'stacked')
        axis([0 61 0 .5])
    end
    %set(gcf,'Units','Normalized','OuterPosition',[0 0 .5 1])
    %saveas(gcf,['/Users/ugne/Dropbox/mouseMASTER/grooms/day_' num2str(d) '.tiff']);
end


for con = 1:8
    for d = 1:4
gt = gcon{con}(:,:,d);
gtm = median(gt,'omitnan');
gsummary(con,:,d) = gtm;
    end
end

for con = 1:8
    for d = 1:4
gt = gcon{con}(:,:,d);
gtm = mean(gt,'omitnan');
gsummarym(con,:,d) = gtm;
    end
end

%% GROOMING COMPARISON HERE
for con = 1:8
    subplot(3,3,con); 
bar(squeeze(gsummary(con,1:5,:))','stacked')
axis([0 5 0 .25]);
end

   figure(2);      
   for con = 1:8
    subplot(3,3,con); 
bar(squeeze(gsummarym(con,1:5,:))','stacked')
axis([0 5 0 .25]);
   end
   
   
cmapgroom = [128 177 211; 251 128 114; 255 255 179; 253 180 98; 141 211 199]./256;
 
gs1 = squeeze(gsummarym(1,1:5,:))';
gs2 = squeeze(gsummarym(2,1:5,:))';
gs3 = squeeze(gsummarym(3,1:5,:))';
gsC2 = [gs1; gs2; gs3];
bx = [1 2 3 4 6 7 8 9 11 12 13 14];
figure(1);
b1 = bar(bx,gsC2,'stacked');
for i = 1:5
    b1(i).FaceColor = 'flat';
    b1(i).CData = cmapgroom(i,:);
end
axis([0 15 0 .25]);
saveas(gcf,['/Users/ugne/Dropbox/mouseMaster/groom_C2.pdf'])
axis off
saveas(gcf,['/Users/ugne/Dropbox/mouseMASTER/groom_C2_NA.pdf'])
close(gcf)



gs1 = squeeze(gsummarym(4,1:5,:))';
gs2 = squeeze(gsummarym(5,1:5,:))';
gs3 = squeeze(gsummarym(6,1:5,:))';
gsL7 = [gs1; gs2; gs3];
bx = [1 2 3 4 6 7 8 9 11 12 13 14];
figure(2);
b2 = bar(bx,gsL7,'stacked');
for i = 1:5
    b2(i).FaceColor = 'flat';
    b2(i).CData = cmapgroom(i,:);
end
axis([0 15 0 .25]);
saveas(gcf,['/Users/ugne/Dropbox/mouseMaster/groom_L7.pdf'])
axis off
saveas(gcf,['/Users/ugne/Dropbox/mouseMASTER/groom_L7_NA.pdf'])
close(gcf)

gs1 = squeeze(gsummarym(7,1:5,:))';
gs2 = squeeze(gsummarym(8,1:5,:))';
gsWT = [gs1; gs2;];
bx = [1 2 3 4 6 7 8 9];
figure(3);
b3 = bar(bx,gsWT,'stacked');
for i = 1:5
    b3(i).FaceColor = 'flat';
    b3(i).CData = cmapgroom(i,:);
end
axis([0 10 0 .25]);
saveas(gcf,['/Users/ugne/Dropbox/mouseMaster/groom_WT.pdf'])
axis off
saveas(gcf,['/Users/ugne/Dropbox/mouseMASTER/groom_WT_NA.pdf'])
close(gcf)   
