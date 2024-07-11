function [m, mse, n] = calcGroupMean(v, g, groups, data_type)

if nargin<3
    groups = unique(g);
end

if nargin<4 || isempty(data_type)
    data_type = 'continuous';
end

if ~iscell(groups)
    m = arrayfun(@(s) mean(v(g==s)), groups);
    n = arrayfun(@(s) sum(g==s), groups);
    switch data_type
        case 'binary'
            mse = arrayfun(@(s) sqrt(m(groups==s).*(1-m(groups==s))./sum(g==s)), groups);
        case 'continuous'
            mse = arrayfun(@(s) std(v(g==s))/sqrt(sum(g==s)), groups);
    end
else
    m = cellfun(@(s) mean(v(ismember(g,s))), groups);
    n = cellfun(@(s) sum(ismember(g,s)), groups);
    switch data_type
        case 'binary'
            findcell = @(c,pat) cellfun(@(x,s) isequal(x,s), c, repmat({pat},size(c)));
            mse = cellfun(@(s) sqrt(m(findcell(groups,s)).*(1-m(findcell(groups,s)))./sum(ismember(g,s))), groups);
        case 'continuous'
            mse = cellfun(@(s) std(v(ismember(g,s)))/sqrt(sum(ismember(g,s))), groups);
    end
end

end