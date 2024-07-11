function [fh, fh_idv, stat] = run_show_dthreshold_logistic(cond_switch, coh, resp, targ_cor, stimdur, subj, opt)

opt.plot_idv_fh = false;
opt.binning_method = 'quantile';

nsubj = length(opt.subj_list);
for s = 1:nsubj
    I = strcmp(subj, opt.subj_list{s});
    opt.subj_idx = s;
    [stat{s}, fh_idv{s}] = show_dthreshold(cond_switch(I), coh(I), resp(I), targ_cor(I), stimdur(I), opt);
end
fh = show_dthreshold_average(stat, opt);

end


function fh = show_dthreshold_average(stat, opt)

switch opt.stimdur_or_csi
    case 'stimdur'  
        opt.xlabel = 'Stimulus duration (ms)';
    case 'csi'
        opt.xlabel = 'CSI (ms)';
end

nsubj = length(stat);
dur_bin_center = stat{1}.dur_bin_center;

%% data averaging
thres1 = cellfun(@(x) x.thres1, stat, 'uni', 0); thres1 = cell2mat(thres1);
thres1_mean = mean(thres1, 2); thres1_sem = std(thres1, 0, 2) / sqrt(nsubj);
thres2 = cellfun(@(x) x.thres2, stat, 'uni', 0); thres2 = cell2mat(thres2);
thres2_mean = mean(thres2, 2); thres2_sem = std(thres2, 0, 2) / sqrt(nsubj);
dthres = cellfun(@(x) x.dthres, stat, 'uni', 0); dthres = cell2mat(dthres);
dthres_mean = mean(dthres, 2); dthres_sem = std(dthres, 0, 2) / sqrt(nsubj);

%% linear regression line averaging
fdthres = cellfun(@(x) x.fdthres, stat, 'uni', 0); fdthres = cell2mat(fdthres);
fdthres_mean = mean(fdthres, 2); fdthres_sem = std(fdthres, 0, 2) / sqrt(nsubj);

%% plot
fh = figure('color','w','Position',[100 100 350 150]);

% threshold of switch and non-switch trials
subplot(1,2,1)
hold on
plot(dur_bin_center, thres1_mean, 'marker', '.', 'color', 'k', 'markersize', 7);
% cerrorbar(dur_bin_center, thres1_mean, thres1_sem, 'Color', 'k');
plot(dur_bin_center, thres2_mean, 'marker', '.', 'color', 'r', 'markersize', 7);
% cerrorbar(dur_bin_center, thres2_mean, thres2_sem, 'Color', 'r');
set(gca, 'xscale', 'log', 'xminortick','off', ...
    'yscale', 'log', 'yminortick','off');
if strcmp(opt.stimdur_or_csi,'csi')
    xlim([500 1200])
    ylim([15 60])
elseif strcmp(opt.stimdur_or_csi,'stimdur')
    xlim([300 700])
    ylim([15 60])
end
xlabel(opt.xlabel)
ylabel('Threshold (% Morph)')
if strcmp(opt.stimdur_or_csi,'csi')
    xticks(550:150:1150)
    xticklabels({'550','','850','','1150'})
end

end


function [stat, fh] = show_dthreshold(cond, coh, resp, targ_cor, stimdur, opt)

coh = abs(coh);
corr = resp==targ_cor;

%% get stim dur bin
switch opt.stimdur_or_csi
    case 'stimdur'
        quantile = {320,427,533,640};     
        opt.xlabel = 'Stimulus duration (ms)';
    case 'csi'
        quantile = split_epoch(stimdur, opt.binning_method, opt.nbin, false);
        opt.xlabel = 'CSI (ms)';
end
dur_bin = nan(length(stimdur),1);
dur_bin_center = nan(length(quantile),1);
for n = 1:length(quantile)
    idx = ismember(stimdur, quantile{n});
    dur_bin(idx) = n;
    dur_bin_center(n) = exp(mean(log(stimdur(dur_bin==n))));
end

%% get threshold from logistic fit
ucoh = unique(coh);
for n = 1:length(dur_bin_center)
    CORR = corr(dur_bin==n);
    COH = coh(dur_bin==n);
    COND  = cond(dur_bin==n);
    % data averaging
    [pc1{n}, pc1se{n}] = calcGroupMean(CORR(COND==1), COH(COND==1), ucoh, 'binary');
    [pc2{n}, pc2se{n}] = calcGroupMean(CORR(COND==2), COH(COND==2), ucoh, 'binary');
    % get threshold from logistic fit
    alpha{n} = glmfit([COND==2, COH, COH.*(COND==2)], CORR, 'binomial', 'link', 'logit'); % logit(P) = a0 + a0s * cond + (a1 + a1s * cond) * coh
    thres1(n) = fminsearch(@(coh) abs(glmval(alpha{n}, [1==2, coh, coh.*(1==2)], 'logit') - 0.816), .2);
    thres2(n) = fminsearch(@(coh) abs(glmval(alpha{n}, [2==2, coh, coh.*(2==2)], 'logit') - 0.816), .2);
end
thres1 = thres1*100;
thres2 = thres2*100;

%% get dthreshold
dthres = thres2 - thres1;
b = regress(dthres(:), [ones(length(dur_bin_center),1) dur_bin_center(:)]);
fdthres = b(1) + b(2) * dur_bin_center;

%% plot threshold
if opt.plot_idv_fh
    fh = figure('color','w','Position',[200+200*opt.subj_idx 100 200 100]);
    
    % threshold
    subplot(1,2,1);
    hold on;
    plot(dur_bin_center, thres1, 'marker', '.', 'color', 'k', 'markersize', 7);
    plot(dur_bin_center, thres2, 'marker', '.', 'color', 'r', 'markersize', 7);
    set(gca, 'xscale', 'log', 'xminortick','off', ...
        'ylim', [min([thres1 thres2])*.8 max([thres1 thres2])*1.2], 'yscale', 'log', 'yminortick','off');
    if strcmp(opt.stimdur_or_csi,'csi')
        xticks([600 700 800 900 1000])
        xticklabels({'600','','800','','1000'})
    end
    axis square;
    
    % dthreshold
    subplot(1,2,2);
    hold on;
    plot(dur_bin_center, dthres, 'marker', '.', 'color', 'k', 'markersize', 7); % datapoint
    set(gca, 'xscale', 'log', 'xminortick','off')
    axis square;
else
    fh = [];
end

%% save
stat = struct('dur_bin_center', dur_bin_center, 'thres1', thres1(:), 'thres2', thres2(:), ...
    'dthres', dthres(:), 'fdthres', fdthres(:));

end