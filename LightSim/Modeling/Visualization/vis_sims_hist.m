
for lesion_i = 1:19%numel(scaled_metadata)
    lmd = scaled_metadata(lesion_i);
    image((lmd.slice+2*lmd.scaled_simslice)/3)
    daspect([1 1 1])
    disp(lmd.name)

end