% CODAANALYSIS: compositional data analysis of fractions of time spent in
% behaviors; log ratio difference between groups and control group

fig_path = '../Figures/';
data_path = '../Data/';
BehaviorLabels = {'IDLE','GROOM',['SLOW' newline 'EXPLORE'],['FAST' newline 'EXPLORE'],'REAR','CLIMB',['AMBLE/' newline 'TURN'],'LOCOMOTION'};

% if fig_path does not exist, create new folder
if ~exist(fig_path,'dir')
    mkdir(fig_path)
end

%% load the data

probs = readtable([data_path 'MMM_8_data_set.csv']);
% delete NaNs (mouse 16, C57Bl-FEMALE on day 2)
probs = probs(~isnan(table2array(probs(:,5))),:);

day = 4;

probs = probs(probs.day == day,:);
probs = unstack(probs,'prob','class');
probs = probs(:,[1 4:end]);

%% male vs female

probs_sel = probs(strcmp(probs.group,'C57bl-MALE') | strcmp(probs.group,'C57bl-FEMALE'),:);
ref_group = 'C57bl-MALE';
groups = unique(probs_sel.group);
n_groups = numel(groups);

% bootstrapping to get estimate of variability
nBoot = 5000;

bootstat = cell(1,n_groups);
for g = 1:n_groups
    mat = table2array(probs_sel(strcmp(probs_sel.group,groups{g}),2:end));
    bootstat{g} = bootstrp(nBoot,@geomean,mat);
end

idx = cellfun(@(x) strcmp(x,ref_group),groups);
groups_of_interest = find(~idx); % all groups except reference group

CIs_lower_bootstrap = cell(1,length(groups_of_interest));
CIs_upper_bootstrap = cell(1,length(groups_of_interest));
means_bootstrap = cell(1,length(groups_of_interest));
for g = 1:numel(groups_of_interest)
    data1 = bootstat{groups_of_interest(g)};
    data2 = bootstat{idx};
    logratio_tmp = bsxfun(@rdivide,data1,data2); % same as data1./data2
    logratio_tmp = log(logratio_tmp);
    means_bootstrap{g} = mean(logratio_tmp);
    CIs_lower_bootstrap{g} = prctile(logratio_tmp,2.5);
    CIs_upper_bootstrap{g} = prctile(logratio_tmp,97.5);
end
means_bootstrap = vertcat(means_bootstrap{:});
CIs_lower_bootstrap = vertcat(CIs_lower_bootstrap{:});
CIs_upper_bootstrap = vertcat(CIs_upper_bootstrap{:});

% plot
f = figure;
f.Units = 'centimeters';
f.Position = [10,10,13,4];
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)+1])

hb = bar(means_bootstrap(:,1:end)',0.6,'EdgeColor','k','FaceColor','w');
hold on
for i = 1:numel(hb)
    xData = hb(i).XData+hb(i).XOffset;
    line([xData;xData],[CIs_lower_bootstrap(i,:);CIs_upper_bootstrap(i,:)],'Color','k');
end
for i = 1:numel(BehaviorLabels)
    text(i,1.1,BehaviorLabels{i},'HorizontalAlignment','center','FontSize',6)
end
ylim([-2 1])
prepfig()
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'XColor','none')
box off
ax = gca;
ax.FontSize = 9;
set(gca,'XColor','none')
ylabel('ln(Female / Male)','Fontsize',9)

% save the figure
print(gcf,[fig_path 'fraction_log_ratio_female_vs_male_day' num2str(day) '.pdf'],'-dpdf','-r0');
close(gcf)

%% Cntnap2

probs_sel = probs(strcmp(probs.group,'Cntnap2-NEG') | strcmp(probs.group,'Cntnap2-HET') | ...
    strcmp(probs.group,'Cntnap2-HOMO'),:);
ref_group = 'Cntnap2-NEG';
groups = unique(probs_sel.group);
n_groups = numel(groups);

% bootstrapping to get estimate of variability
nBoot = 5000;

bootstat = cell(1,n_groups);
for g = 1:n_groups
    mat = table2array(probs_sel(strcmp(probs_sel.group,groups{g}),2:end));
    bootstat{g} = bootstrp(nBoot,@geomean,mat);
end

idx = cellfun(@(x) strcmp(x,ref_group),groups);
groups_of_interest = find(~idx); % all groups except reference group

CIs_lower_bootstrap = cell(1,length(groups_of_interest));
CIs_upper_bootstrap = cell(1,length(groups_of_interest));
means_bootstrap = cell(1,length(groups_of_interest));
for g = 1:numel(groups_of_interest)
    data1 = bootstat{groups_of_interest(g)};
    data2 = bootstat{idx};
    logratio_tmp = bsxfun(@rdivide,data1,data2); % same as data1./data2
    logratio_tmp = log(logratio_tmp);
    means_bootstrap{g} = mean(logratio_tmp);
    CIs_lower_bootstrap{g} = prctile(logratio_tmp,2.5);
    CIs_upper_bootstrap{g} = prctile(logratio_tmp,97.5);
end
means_bootstrap = vertcat(means_bootstrap{:});
CIs_lower_bootstrap = vertcat(CIs_lower_bootstrap{:});
CIs_upper_bootstrap = vertcat(CIs_upper_bootstrap{:});

% plot
f = figure;
f.Units = 'centimeters';
f.Position = [10,10,13,7];
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)+1])

t = tiledlayout(2,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for g = 1:(n_groups-1)
    nexttile
    hb = bar(means_bootstrap(g,1:end)',0.6,'EdgeColor','k','FaceColor','w');
    hold on
    xData = hb.XData+hb.XOffset;
    line([xData;xData],[CIs_lower_bootstrap(g,:);CIs_upper_bootstrap(g,:)],'Color','k');
    if g == 1
        for i = 1:numel(BehaviorLabels)
            text(i,1.1,BehaviorLabels{i},'HorizontalAlignment','center','FontSize',6)
        end
    end
    ylim([-2 1])
    prepfig()
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'XColor','none')
    box off
    ax = gca;
    ax.FontSize = 9;
    set(gca,'XColor','none')
    ylabel(['ln(' groups{groups_of_interest(g)} ' / ' newline ref_group ')'],'Fontsize',9)
end

% save the figure
print(gcf,[fig_path 'fraction_log_ratio_cntnap2_day' num2str(day) '.pdf'],'-dpdf','-r0');
close(gcf)

%% L7-Tsc1

probs_sel = probs(strcmp(probs.group,'L7-Tsc1-HOM') | strcmp(probs.group,'L7-Tsc1-HET') | ...
    strcmp(probs.group,'L7-Tsc1-NEG'),:);
ref_group = 'L7-Tsc1-NEG';
groups = unique(probs_sel.group);
n_groups = numel(groups);

% bootstrapping to get estimate of variability
nBoot = 5000;

bootstat = cell(1,n_groups);
for g = 1:n_groups
    mat = table2array(probs_sel(strcmp(probs_sel.group,groups{g}),2:end));
    bootstat{g} = bootstrp(nBoot,@geomean,mat);
end

idx = cellfun(@(x) strcmp(x,ref_group),groups);
groups_of_interest = find(~idx); % all groups except reference group

CIs_lower_bootstrap = cell(1,length(groups_of_interest));
CIs_upper_bootstrap = cell(1,length(groups_of_interest));
means_bootstrap = cell(1,length(groups_of_interest));
for g = 1:numel(groups_of_interest)
    data1 = bootstat{groups_of_interest(g)};
    data2 = bootstat{idx};
    logratio_tmp = bsxfun(@rdivide,data1,data2); % same as data1./data2
    logratio_tmp = log(logratio_tmp);
    means_bootstrap{g} = mean(logratio_tmp);
    CIs_lower_bootstrap{g} = prctile(logratio_tmp,2.5);
    CIs_upper_bootstrap{g} = prctile(logratio_tmp,97.5);
end
means_bootstrap = vertcat(means_bootstrap{:});
CIs_lower_bootstrap = vertcat(CIs_lower_bootstrap{:});
CIs_upper_bootstrap = vertcat(CIs_upper_bootstrap{:});

% plot
f = figure;
f.Units = 'centimeters';
f.Position = [10,10,13,7];
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)+1])

t = tiledlayout(2,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for g = 1:(n_groups-1)
    nexttile
    hb = bar(means_bootstrap(g,1:end)',0.6,'EdgeColor','k','FaceColor','w');
    hold on
    xData = hb.XData+hb.XOffset;
    line([xData;xData],[CIs_lower_bootstrap(g,:);CIs_upper_bootstrap(g,:)],'Color','k');
    if g == 1
        for i = 1:numel(BehaviorLabels)
            text(i,1.1,BehaviorLabels{i},'HorizontalAlignment','center','FontSize',6)
        end
    end
    ylim([-2 1])
    prepfig()
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'XColor','none')
    box off
    ax = gca;
    ax.FontSize = 9;
    set(gca,'XColor','none')
    ylabel(['ln(' groups{groups_of_interest(g)} ' / ' newline ref_group ')'],'Fontsize',9)
end

% save the figure
print(gcf,[fig_path 'fraction_log_ratio_L7Tsc1_day' num2str(day) '.pdf'],'-dpdf','-r0');
close(gcf)