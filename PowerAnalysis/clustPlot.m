function clustPlot(powerArray,groupings,colors)
% function to plot power vs time with color-coded channel Groupings

channels = 1:32;
alpha = 0.2;

for iChan = channels
    clust = groupings(iChan); 
    if ~isnan(clust)
        plot(powerArray(:,iChan), 'linewidth', 2, 'color', [colors(clust,:) alpha]), hold on
    end
end
end