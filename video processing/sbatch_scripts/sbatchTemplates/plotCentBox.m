set(0,'DefaultFigureVisible','off')
file_dir=[pwd '/box/'];
if ~isdir(file_dir)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', file_dir)
    return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(file_dir, '*box_data.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(file_dir, baseFileName);
    fprintf(1, 'file %s\n', fullFileName);
end

%%
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(file_dir, baseFileName);
    load(fullFileName);
    c = linspace(1,10,length(mouseDataBox.centroid));
h=figure; set(h,'visible','off');
    scatter(mouseDataBox.centroid(:,1),mouseDataBox.centroid(:,2), 2, c , 'filled')
    set(gcf,'renderer','Painters') % rendering for eps-files
    set(gca,'PlotBoxAspectRatio',[1 1 1]) % Aspect Rario
    set(gca,'Color','w')
    set(gcf,'Color','w') % match figure background
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
        'PaperUnits', 'Inches', 'PaperSize', [5, 5]) % size and location in inches
    box off; axis off
    title(strrep(baseFileName(1:end-4), '_', ' '), 'FontSize', 14);
    print(gcf,[file_dir, baseFileName(1:end-4),'.png'],'-dpng','-r300');
    close all
end
