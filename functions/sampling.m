function [edge, csi_bin_median] = sampling(x, opt)

if strcmp(opt.sample,'log')
    % cumulative distribution of csi
    csi_sort = sort(x); % in ascending order
    csi_iso_samp = linspace(min(x),max(x),100)'; % isometric sampling of csi
    csi_cum_distribution = zeros(100,1); % cumulative distribution of csi
    for k = 1:length(csi_cum_distribution)
        I = csi_sort<=csi_iso_samp(k);
        csi_cum_distribution(k) = sum(I) / length(I);
    end
    
    % iso-area sampling (a method taking into account substitution)
    p = linspace(0,1,opt.nbin+1)'; % prob
    edge = nan(opt.nbin+1, 1); % edges of sampling bin
    for k = 1:length(p)
        [~, idx] = min(abs(p(k)-csi_cum_distribution));
        edge(k) = csi_iso_samp(idx);
    end
    csi_bin_median = (edge(1:end-1) + edge(2:end)) / 2; % median of each bin
end

if strcmp(opt.sample,'linear')
    edge = linspace(min(x), max(x), opt.nbin+1)';
    csi_bin_median = (edge(1:end-1) + edge(2:end)) / 2; % median of each bin
end

end