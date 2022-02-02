function vis_bounds(bounds)
    X = 0:.053:.053*(size(bounds,1)-1);
    Y = 0:.053:.053*(size(bounds,2)-1);
    Z = 0:.45:.45*(size(bounds,3)-1);
    
    figure;
    isosurface(Y, X, Z, bounds);
    daspect([1 1 1]);
end