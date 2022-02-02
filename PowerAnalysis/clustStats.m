function [p,table] = clustStats(clust1,clust2,clust3)

if isempty(clust1) && isempty(clust2)
    warning('No multiple comparisons possible')
    p = NaN; table = cell(3,6);
    return
elseif isempty(clust1) && isempty(clust3)
    warning('No multiple comparisons possible')
    p = NaN; table = cell(3,6);
    return
elseif isempty(clust2) && isempty(clust3)
    warning('No multiple comparisons possible')
    p = NaN; table = cell(3,6);
    return
end

allClusts = [clust1 clust2 clust3]; % combine clusters to 1 array

Label = cell(size(allClusts));
for i = 1:length(clust1)
    Label{1,i} = 'Low';
end
for i = (length(clust1)+1):(length(clust1)+length(clust2))
    Label{1,i} = 'Medium';
end
for i = (length(clust2)+length(clust1)+1):length(allClusts)
    Label{1,i} = 'High';
end

% one-way anova test
[p,anovatab,s] =  anova1(allClusts,Label,'off');

% multiple comparison between groups
[c,m,h,nms] = multcompare(s,'display','off');

% table of comparisons
table = [nms(c(:,1)), nms(c(:,2)), num2cell(c(:,3:6))];

if isempty(clust1) || isempty(clust2) || isempty(clust3)
    for i = 1:6
        table{2,i} = NaN; 
        table{3,i} = NaN;
    end
    
end

end