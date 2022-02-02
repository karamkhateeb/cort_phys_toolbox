function simslice = rescale_sims(simslice)
    [newsimdepth, newsimwidth] = size(simslice);
    
    %make sure they reach to 800 and 500
    if newsimdepth < 500
        depth_add = 500 - newsimdepth;
        simslice = cat(1, simslice, zeros(depth_add, newsimwidth));
    else
        simslice = simslice(1:500, :);
    end
    if newsimwidth < 800
        width_add = 800 - newsimwidth;
        simslice = cat(2, zeros(500, floor(width_add/2)), simslice, zeros(500, ceil(width_add/2)));
    else
        simslice = simslice(:, 1:800);
    end
end
