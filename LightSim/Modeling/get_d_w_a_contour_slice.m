function [maxdepth, avgwidth, area, midwidth] = get_d_w_a_contour_slice(contour)
    %% Get metrics for a contour
    model_c_mask = imrotate(256*poly2mask(contour.xdata, contour.ydata, 500, 800), 180);
    [maxdepth, avgwidth, area, midwidth] = get_lesion_maxdepths_and_avgwidths(model_c_mask);
end