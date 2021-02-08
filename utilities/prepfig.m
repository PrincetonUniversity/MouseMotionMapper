function prepfig()
% set some graphical attributes of the current axis

set(gcf,'renderer','Painters') % rendering for eps-files
set(gca,'TickDir','out') % draw the tick marks on the outside
set(gca,'PlotBoxAspectRatio',[1 1 1]) % Aspect Rario
set(gca,'Color','w')
set(gcf,'Color','w') % match figure background
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
    'PaperUnits', 'Inches', 'PaperSize', [5, 5]) % size and location in inches
box off
set(gca,'FontSize',12) % Creates an axes and sets its FontSize to 18
set(gca, 'ycolor', 'k') % black axis
set(gca, 'xcolor', 'k')