% first load stuff in memory
full_metadata = scaled_metadata;
c_d_w_a = zeros(numel(full_metadata), 4);
for lesion_i = 1:numel(full_metadata)
    lmd = full_metadata(lesion_i);
    
%     switch lmd.aperture
%         case 0.5
%             model_C_data = model_C_data_p5mm;
%         case 1
%             model_C_data = model_C_data_1mm;
%         case 1.5
%             model_C_data = model_C_data_1p5mm;
%         case 2
%             model_C_data = model_C_data_2mm;
%     end
        
%     cntr = get_model_contour(lmd.aperture, lmd.intensity, 7.9e-05, model_C_data); %Here
%     cntr = get_model_contour(lmd.aperture, lmd.intensity, 0.00013301, model_C_data); %Here
%     cntr = get_model_contour(lmd.aperture, lmd.intensity, 0.00019862, model_C_data); %Here
    cntr = lmd.simslice;

    [depth, width, area, midwidth] = get_lesion_maxdepths_and_avgwidths(cntr);
%     [depth, width, area, midwidth] = get_d_w_a_contour_slice(cntr);
    c_d_w_a(lesion_i, :) = [depth, width, area, midwidth];
end