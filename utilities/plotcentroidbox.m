close all, clear all

%% Check to make sure that folder actually exists.  Warn user if it doesn't.
File_dir='/tigress/jverpeut/FMB_DREADD4/box/';
if ~isdir(File_dir)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', File_dir);
    uiwait(warndlg(errorMessage));
    return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(File_dir, '*box_data.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(File_dir, baseFileName);
    fprintf(1, 'file %s\n', fullFileName);
end

%%
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(File_dir, baseFileName);
    load(fullFileName);
    c = linspace(1,10,length(mouseDataBox.centroid));
    figure
    scatter(mouseDataBox.centroid(:,1),mouseDataBox.centroid(:,2), 2, c , 'filled')
    set(gcf,'renderer','Painters') % rendering for eps-files
    set(gca,'PlotBoxAspectRatio',[1 1 1]) % Aspect Rario
    set(gca,'Color','w')
    set(gcf,'Color','w') % match figure background
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
        'PaperUnits', 'Inches', 'PaperSize', [5, 5]) % size and location in inches
    box off; axis off
    title(strrep(baseFileName(1:end-4), '_', ' '), 'FontSize', 14);
    print(gcf,[File_dir, baseFileName(1:end-4),'.png'],'-dpng','-r300');
    close all
end