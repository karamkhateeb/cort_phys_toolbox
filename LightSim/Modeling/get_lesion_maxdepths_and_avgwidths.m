function [maxdepth, avgwidth, area, midwidth] = get_lesion_maxdepths_and_avgwidths(lesion)
    if lesion == 0
        maxdepth = 0;
        avgwidth = 0;
        area = 0;
        midwidth = 0;
    else
        [rowIdcs, colIdcs] = find(lesion~=0);
        maxdepth = max(rowIdcs);
        area = nnz(lesion);
        avgwidth = area / maxdepth;
    %     width = max(colIdcs) - min(colIdcs);
        midline = floor(maxdepth/2);
        midline_vals = colIdcs(rowIdcs == midline);
        midwidth = max(midline_vals) - min(midline_vals);
    end
end