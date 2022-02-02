function clustArray(varargin)
% Use this function to generate a heat map of the ECoG array groupings

narginchk(2,3);
numArgs = nargin;
arrayVals = varargin{1};
map = varargin{2};
if numArgs == 3
    ax = varargin{3};
end

% X and Y locations for array
x = [4 4 3 2 2 2 1 1 1 1 2 2 3 3 3 4 4 4 5 6 6 7 8 7 7 7 6 6 6 5 5 5]*0.75;
y = [3 1 2 1 5 3 4 2 8 6 9 7 6 4 8 5 7 9 8 9 7 8 5 6 2 4 3 5 1 6 4 2]*0.75;
[xq, yq] = meshgrid(0.75:0.05:8*0.75,0.75:0.05:9*0.75);

method = 'nearest';

a = griddata(x, y, arrayVals, xq, yq, method);
pcolor(xq, yq, a), shading flat
if numArgs == 3, colormap(ax,map); else, colormap(map); end
hold on; scatter(x,y,18,'k','filled');
xlim([0.5, 8*0.75+.25]); ylim([0.5, 9*0.75+.25]); box off
set(gca,'ytick',[],'xtick',[], 'xticklabel', [], 'Yticklabel', [], 'DataAspectRatio', [1 1 1]);

end