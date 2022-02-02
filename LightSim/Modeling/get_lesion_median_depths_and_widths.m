function [meddepth, medwidth, area, midwidth] = get_lesion_median_depths_and_widths(lesion)
    [rowIdcs, colIdcs] = find(lesion~=0);
    depths = zeros(1,numel(colIdcs));
    for i = 1:numel(colIdcs)
        col = colIdcs(i);
        [rows,~] = find(lesion(:,col)~=0);
        depths(i) = max(rows);
    end
    meddepth = median(depths);
    maxdepth = max(rowIdcs);
    area = nnz(lesion);
    avgwidth2 = area / maxdepth;
    for j = 1:numel(rowIdcs)
        row = rowIdcs(j);
        [~,cols] = find(lesion(row,:)~=0);
        widths(row) = max(cols) - min(cols);
    end
    medwidth = median(widths);
%     width = max(colIdcs) - min(colIdcs);
    midline = floor(maxdepth/2);
    midline_vals = colIdcs(rowIdcs == midline);
    midwidth = max(midline_vals) - min(midline_vals);
end