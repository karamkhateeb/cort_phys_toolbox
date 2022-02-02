% % Find best threshold to explain all lesions
% % generate an average IOU over slices for each threshold
get_all_slices;
model_suffix = '_gua_yus';
load_all_contours;
load_lesion_metadata;
disp('done loading')

%%
mean_IOUs = [];
for threshold = [1.15] * 1e-4
    disp(threshold)
    IOUs = [];
    for les_i = 1:numel(lesion_metadata)
        lmd = lesion_metadata(les_i);
        cntr = get_model_contour(lmd.aperture, lmd.intensity, threshold);
        IOU = get_IOU_contour_slice(cntr, lmd.slice);
        IOUs = [IOUs IOU];
    end

    mean_IOUs = [mean_IOUs, mean(IOUs)];    
end