function show_IOU(contour, slice)
    model_c_mask = imrotate(256*poly2mask(contour.xdata, contour.ydata, 500, 800), 180);
    
    % show whassup
    image((slice+2*model_c_mask)/3)
    daspect([1 1 1])
end