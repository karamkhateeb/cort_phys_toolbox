% full_metadata = lesion_metadata;
for lesion_i = 1:numel(full_metadata)
    lmd = full_metadata(lesion_i);
    disp(lmd.name)
    
    switch lmd.aperture
        case 0.5
            model_C_data = model_C_data_p5mm;
        case 1
            model_C_data = model_C_data_1mm;
        case 1.5
            model_C_data = model_C_data_1p5mm;
        case 2
            model_C_data = model_C_data_2mm;
    end
        
    cntr = get_model_contour(lmd.aperture, lmd.intensity, 0.00019862, model_C_data);
    model_c_mask = imrotate(256*poly2mask(cntr.xdata, cntr.ydata, 500, 800), 180);
    full_metadata(lesion_i).simslice = model_c_mask;
end