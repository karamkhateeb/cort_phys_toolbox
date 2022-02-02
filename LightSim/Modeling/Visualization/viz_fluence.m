%clear

% data = loadmc2("/home/julien/Documents/MCXStudio/Output/mcxsessions/Cyl_iso/Cyl_iso.mc2", [200, 200, 200]);
% data = loadmc2("C:\Users\julienb\Documents\MCXStudio\Output\mcxsessions\Cyl\Cyl.mc2", [200, 200, 200]);
% data = loadmc2("C:\Users\julienb\Documents\MCXStudio\Output\mcxsessions\560nm_1mm_w_pdms\560nm_1mm_w_pdms.mc2", [80, 80, 57]);
% data = loadmc2("C:\Users\julienb\Documents\MCXStudio\Output\mcxsessions\560nm_1mm\560nm_1mm.mc2", [80, 80, 51]);
% data = loadmc2("old_same/560nm_1mm_w_pdms/560nm_1mm_w_pdms.mc2", [80, 80, 57]);
data = loadmc2("/home/julien/lab/MCXStudio/Output/imported/560nm_1p5mm_w_pdms_e7.mc2", [800, 800, 563]);

log_data = log(data(:, :, 62:end));

% smol_data = log_data(61:140, 61:140, 1:80);
% [X, Y, Z]  = meshgrid(1:100);
% slice(X, Y, Z, log_data, 21, 21, 1)
% heatmap(imrotate(squeeze(smol_data(21, :, :)), 90), 'XLabel', 'X (mm)', 'Title','log(intensity) at central slice: 560nm', ...
%     'YLabel', 'Z (mm)', 'GridVisible', 'off', 'XData', [-4.0:.1:3.9], 'YData', flip([0:.1:7.9]));


%% For 3d contour surface plots
figure
for intens = 3+ceil(min(min(min(log_data)))):.5:floor(max(max(max(log_data))));
    isosurface(log_data, intens);
end


%% For a contour plot
figure
% K = ones(3);
% smooth_log_data = conv2(squeeze(log_data(40, :, :)), K,'same');
smooth_log_data = squeeze(log_data(40, :, :));
[M, c] = contour(imrotate(smooth_log_data, 90), 18);
c.LineWidth = 2;

set(gca,'XTick',0:10:80,'XTickLabel',-4:1:4)
pbaspect([9 5 1]);
set(gca,'YTick',0:10:50,'YTickLabel', flip(0:1:5));

%% For a heatmap plot of data
figure;
h = heatmap(imrotate(squeeze(log_data(:, 40, :)), -90),   ...
     'GridVisible', 'off', 'XData', [-4:.1:3.9], 'YData', [0:.1:4.9]);
% set(gca, 'XTick', [0:10:80]); 
% set(gca, 'XTickLabel', [-40:10:40]);

h.Position(3:4) = min(h.Position(3:4))*[1,1]; 
% 
% % %%
% 
% ylabel('Z (mm)')
% xlabel('X (mm)')
% set(gco, 'Title','log(intensity) at central slice: 560nm')

colormap default
% colorbar
