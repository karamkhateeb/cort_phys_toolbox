function mapArray(varargin)
% Use this function to generate a heat map of the ECoG array power values

resolFactor = 20; % this is to make the diagnoal lines look smoother in the figure

% X and Y locations for array
x = [4 4 3 2 2 2 1 1 1 1 2 2 3 3 3 4 4 4 5 6 6 7 8 7 7 7 6 6 6 5 5 5]*0.75*resolFactor;
y = [3 1 2 1 5 3 4 2 8 6 9 7 6 4 8 5 7 9 8 9 7 8 5 6 2 4 3 5 1 6 4 2]*0.75*resolFactor;
[xq, yq] = meshgrid(0.75:0.05:8*0.75,0.75:0.05:9*0.75);
xq = xq.*resolFactor; yq = yq.*resolFactor;

narginchk(1,3);

numArgs = nargin;

arrayVals = varargin{1};

method = 'nearest'; % 'Linear', 'natural', or 'cubic'

% color map RGB values
map = [linspace(248,64,64) linspace(64,8,64);
        linspace(117,78,64) linspace(78,178,64);
        linspace(117,77,64) linspace(77,227,64)]' / 255;
    
if numArgs == 2 || numArgs == 3    
    if strcmp(varargin{2},'neg')
        map = map(1:64,:);
    elseif strcmp(varargin{2},'pos')
        map = map(64:end,:);
    end
end

if numArgs == 3
    ax = varargin{3};
end

a = griddata(x, y, arrayVals, xq, yq, method);
pcolor(xq, yq, a), shading interp,
if numArgs == 3
    colormap(ax,map);
else
    colormap(map)
end
hold on; scatter(x,y,18,'k','filled');
xlim([0.5, 8*0.75+.25].*resolFactor); ylim([0.5, 9*0.75+.25].*resolFactor); box off
set(gca,'ytick',[],'xtick',[], 'xticklabel', [], 'Yticklabel', [], 'DataAspectRatio', [1 1 1]);

end