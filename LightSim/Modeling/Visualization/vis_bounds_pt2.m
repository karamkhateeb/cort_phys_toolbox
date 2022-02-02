function vis_bounds_pt2(bounds)
    % lesion 1 bounds is p2(200:600,200:1100,27:33)

    X = 0:.025:.025*(size(bounds,1)-1);
    Y = 0:.025:.025*(size(bounds,2)-1);
    Z = 0:.5:.5*(size(bounds,3)-1);
    
    figure;
    isosurface(Y, X, Z, bounds);
    daspect([1 1 1]);
end