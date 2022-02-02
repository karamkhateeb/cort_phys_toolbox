function xpos = getGroupPos(groupedData)
ngroups = size(groupedData, 1); nbars = size(groupedData, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
xpos = nan(1,nbars);
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    xpos = [xpos x];
end
% hold off
end