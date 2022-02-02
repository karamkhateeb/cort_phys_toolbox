scaled_metadata = full_metadata;
for lesion_i = 1:numel(full_metadata)
    lmd = full_metadata(lesion_i);
    simslice = lmd.simslice;
    simslice = imresize(simslice, [500*depth_scaling(lesion_i), 800*width_scaling(lesion_i)], 'nearest');
    simslice = rescale_sims(simslice);
    scaled_metadata(lesion_i).scaled_simslice = simslice;
end