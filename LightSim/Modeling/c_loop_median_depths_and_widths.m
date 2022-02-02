% first load stuff in memory
c_d_w_a = zeros(numel(full_metadata), 4);
for lesion_i = 1:numel(full_metadata)
    lmd = full_metadata(lesion_i);
%     cntr = get_model_contour(lmd.aperture, lmd.intensity, 7.9e-05);
    cntr = lmd.simslice;
    [meddepth, medwidth, area, midwidth] = get_lesion_median_depths_and_widths(cntr);
%     [depth, width, area, midwidth] = get_d_w_a_contour_slice(cntr);
    c_d_w_a(lesion_i, :) = [meddepth, medwidth, area, midwidth];
end