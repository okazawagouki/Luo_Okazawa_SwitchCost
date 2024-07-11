function [quantile, legstr, ID, epoch_mean] = split_epoch(dur, method, nbin, verbose)

switch method
    case 'cumulate'
        [udur, ~, uID] = unique(dur);
        Ndata = histc(uID, 1:max(uID));
        quantile = cell(nbin,1);
        N = zeros(nbin,1);
        idx = 1;
        for m=1:length(Ndata)
            if isempty(quantile{idx})
                quantile{idx} = udur(m);
                N(idx) = Ndata(m);
            elseif N(idx) > sum(N)/nbin
                idx = idx + 1;
                quantile{idx} = udur(m);
                N(idx) = Ndata(m);
            else
                quantile{idx} = [quantile{idx}; udur(m)];
                N(idx) = N(idx) + Ndata(m);
            end
        end
        N(cellfun(@(x)isempty(x), quantile)) = [];
        quantile(cellfun(@(x)isempty(x), quantile)) = [];
        
        ID = nan(length(dur),1);
        for n=1:length(quantile)
            ID(ismember(dur, quantile{n})) = n;
        end
        
        if verbose
            figure('color', 'w', 'position', [100 100 300 300]);
            bar(1:length(N), N, 'hist');
            axis square;
        end

    case 'quantile'
        sdur = sort(dur);
        udur = unique(dur);
        qborder = nan(nbin-1,1);
        for n=1:length(qborder)
            qborder(n) = sdur(round(length(sdur)/nbin * n));
        end
        qborder = [min(dur); qborder(:); max(dur)+1];
        quantile = cell(nbin,1);
        ID = nan(length(dur),1);
        for n=1:nbin
            quantile{n} = udur(udur >= qborder(n) & udur < qborder(n+1));
            ID(dur >= qborder(n) & dur < qborder(n+1)) = n;
        end
end

epoch_mean = nan(max(ID),1);
for n=1:max(ID)
    epoch_mean(n) = mean(dur(ID==n));
end



legstr = cell(length(quantile),1);
for n=1:length(quantile)
    if length(quantile{n}) == 1
        legstr{n} = num2str(quantile{n});
    else
        legstr{n} = sprintf('%d-%d', min(quantile{n}), max(quantile{n}));
    end
end


