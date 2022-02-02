function [hl,pl,sortArr] = SortChannelsNeg(varargin)
%SortChannels sorts channels into groups depending on whether or not
%they're decreasing in power.
%   INPUTS:   
%   post is the time-binned post-PT power (nBins x number of Channels)
%   base is the time-binned baseline power (nBins x number of Channels)
%   alpha is the alpha value [0 1] for specifying significance level
%   (optional). Default value is 0.05 if not specified.
%
%   OUTPUTS:
%   hl is 1 if null hypothesis is rejected for left-tailed ttest, else 0
%   pl are the corresponding pvalues for each ttest
%   sortArr is a 1xN array where N is the number of channels describing the
%   identity of each channel

narginchk(2,3)

if nargin == 2
    alpha = 0.05;
elseif nargin == 3
    alpha = varargin{3};
end

post = varargin{1};
base = varargin{2};


% left-tailed test
hl = NaN(size(base,2),1); pl = NaN(size(base,2),1);
for iChan = 1:size(base,2)
    [hl(iChan), pl(iChan)] = ttest(post(:,iChan),base(:,iChan),'Tail','left','alpha',alpha);
end


% make clustered channel array
sortArr = NaN(size(base,2),1);
sortArr(find(hl == 1)) = 1; % if channel is negative, ID = 1; else, ID = 2
sortArr(find(hl == 0)) = 2;

end

