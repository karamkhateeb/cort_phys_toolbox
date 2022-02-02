% first load stuff in memory
full_metadata = scaled_metadata;
d_w_a = zeros(numel(full_metadata), 4);
for lesion_i = 1:numel(full_metadata)
    lmd = full_metadata(lesion_i);
    lesion = lmd.slice;
    if lesion == 0
        [depth, width, area, midwidth] = deal(0, 0, 0, 0);
    else
        [depth, width, area, midwidth] = get_lesion_maxdepths_and_avgwidths(lesion);
    end
    d_w_a(lesion_i, :) = [depth, width, area, midwidth];
end