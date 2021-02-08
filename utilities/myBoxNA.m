function [outinfo,tt] = myBoxNA(C,mids,w,cvcon,addPoints)

% C: cell array where each cell has an array of measurements
% clr: colors for each array in C
% outinfo: median, 
% find axis sizes - max and min, or set to zero
cc = combineCells(C);
mincc = min(cc); maxcc = 1.1*max(cc);
if mincc<0
    mincc = mincc*1.1;
end
% 
lc = length(C);


%outinfo = zeros(lc,4);

for i = 1:lc
points = C{i};
points = points(~isnan(points));
% whiskers 'min' and 'max' based on inter-quartile range
mp = median(points,'omitnan'); q1 = quantile(points',.25); q3 = quantile(points',.75); 
iqr = q3-q1;
mmin = q1-1.5*iqr;
mmax = q3+1.5*iqr;

if min(points) > mmin
    mmin = min(points);
end
if max(points) < mmax
    mmax = max(points);
end

% points
if addPoints
   np = length(points);
   xl = linspace(mids(i)-w/3,mids(i)+w/3,np);
   scatter(xl,sort(points),5,'filled','MarkerEdgeColor','k','MarkerFaceColor',cvcon(i,:),'LineWidth',1.5);
    
else
% outliers 
under = find(points<mmin);
over = find(points>mmax);
out = [points(under); points(over)];
if ~isempty(out)
    scatter(repmat(mids(i),[1 length(out)]),out,10,'filled','MarkerFaceColor',cvcon(i,:));
end

end

outinfo(i,1) = mp; outinfo(i,2) = q1; outinfo(i,3) = q3; outinfo(i,4) = iqr;
% quantile box (25:75)
corners = [(mids(i)-w/2) (mids(i)+w/2) q1 q3];
fourpt = [corners(1) corners(3) corners(2)-corners(1) corners(4)-corners(3)];
rectangle('Position',fourpt,'EdgeColor',cvcon(i,:),'LineWidth',2);
hold on;


% median line
line([mids(i)-w/1.5 mids(i)+w/1.5],[mp mp],'Color',cvcon(i,:),'LineWidth',3);

line([mids(i) mids(i)],[corners(3) mmin],'LineStyle','--','Color','k','LineWidth',2);
line([mids(i) mids(i)],[corners(4) mmax],'LineStyle','--','Color','k','LineWidth',2);

line([mids(i)-w/4 mids(i)+w/4],[mmin mmin],'Color',cvcon(i,:),'LineWidth',3);
line([mids(i)-w/4 mids(i)+w/4],[mmax mmax],'Color',cvcon(i,:),'LineWidth',3);

end

% statistical tests
tt = zeros(lc,lc);
for i = 1:lc
    for j = 1:lc
        [temph,tempp] = ttest2(C{i},C{j});
        tt(i,j) = tempp;
    end
end


bincoeff = nchoosek(lc,2);
tt = tt.*bincoeff;


% for i = 1:startn
% line([mids(allidx(i)) mids(allidy(i))],[ys(i) ys(i)],'LineWidth',allstar(i),'Color','k')
% scatter(mids(allidx(i)),ys(i),50,cvcon(allidx(i),:),'filled');
% scatter(mids(allidy(i)),ys(i),50,cvcon(allidy(i),:),'filled');
% 
% %mdpt = (mids(allidx(i))+mids(allidy(i)))/2;
% 
% % if allstar(i)==1
% %     scatter(mdpt,ys(i),100,'*','k')
% % elseif allstar(i)==2
% %     scatter(mdpt,ys(i),100,'*','k')
% % elseif allstart(i)==3
% %     scatter(mdpt,ys(i),100,'*','k')
% % end
% end
if mincc>0 
    mincc = 0;
else
    mincc = mincc-.05*(maxcc-mincc);
end
axis([0 lc+1 mincc 1.25*maxcc])








