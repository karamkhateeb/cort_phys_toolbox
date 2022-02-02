for i = 1:size(pt3_L6, 3)
    disp(i)
    disp(sum(sum(pt3_L6(:,:,i))))
    if sum(sum(pt3_L6(:,:,i))) > 660736
        figure
        image(256*pt3_L6(:,:,i))
    end
end