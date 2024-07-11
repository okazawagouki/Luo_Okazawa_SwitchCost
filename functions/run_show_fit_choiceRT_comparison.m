function [fh, fh_idv, stat] = run_show_fit_choiceRT_comparison(coh_ns, resp_ns, fit_pchoice1_ns, rt_ns, fit_rt_ns, subj_ns, ...
    coh_s, resp_s, fit_pchoice1_s, rt_s, fit_rt_s, subj_s, ...
    opt)

opt.plot_idv_fh = false;
opt.RT_stat = 'mean'; % mean or median

nsubj = length(opt.subj_list);
for s = 1:nsubj
    I1 = strcmp(subj_ns, opt.subj_list{s});
    I2 = strcmp(subj_s,  opt.subj_list{s});
    opt.idv_subj = opt.subj_list{s}; opt.idv_param_final = opt.param_final{s};
    [fh_idv{s}, stat{s}] = show_fit_choiceRT(coh_ns(I1), resp_ns(I1), fit_pchoice1_ns(I1), rt_ns(I1), fit_rt_ns(I1), ...
        coh_s(I2),  resp_s(I2),  fit_pchoice1_s(I2),  rt_s(I2),  fit_rt_s(I2), ...
        opt);
end
fh = show_fit_choiceRT_average(stat, opt);

end


function fh = show_fit_choiceRT_average(stat, opt)

nsubj = length(stat);
ucoh = stat{1}.ucoh;

% averaging
p_ns = cellfun(@(x) x.p_ns, stat, 'uni', 0);           p_ns = cell2mat(p_ns);           p_mean_ns = mean(p_ns, 2);           p_sem_ns = std(p_ns, 0, 2) / sqrt(nsubj);
fit_p_ns = cellfun(@(x) x.fit_p_ns, stat, 'uni', 0);   fit_p_ns = cell2mat(fit_p_ns);   fit_p_mean_ns = mean(fit_p_ns, 2);   fit_p_sem_ns = std(fit_p_ns, 0, 2) / sqrt(nsubj);
rtm_ns = cellfun(@(x) x.rtm_ns, stat, 'uni', 0);       rtm_ns = cell2mat(rtm_ns);       rtm_mean_ns = mean(rtm_ns, 2);       rtm_sem_ns = std(rtm_ns, 0, 2) / sqrt(nsubj);
fit_rt_ns = cellfun(@(x) x.fit_rt_ns, stat, 'uni', 0); fit_rt_ns = cell2mat(fit_rt_ns); fit_rt_mean_ns = mean(fit_rt_ns, 2); fit_rt_sem_ns = std(fit_rt_ns, 0, 2) / sqrt(nsubj);

p_s = cellfun(@(x) x.p_s, stat, 'uni', 0);           p_s = cell2mat(p_s);           p_mean_s = mean(p_s, 2);           p_sem_s = std(p_s, 0, 2) / sqrt(nsubj);
fit_p_s = cellfun(@(x) x.fit_p_s, stat, 'uni', 0);   fit_p_s = cell2mat(fit_p_s);   fit_p_mean_s = mean(fit_p_s, 2);   fit_p_sem_s = std(fit_p_s, 0, 2) / sqrt(nsubj);
rtm_s = cellfun(@(x) x.rtm_s, stat, 'uni', 0);       rtm_s = cell2mat(rtm_s);       rtm_mean_s = mean(rtm_s, 2);       rtm_sem_s = std(rtm_s, 0, 2) / sqrt(nsubj);
fit_rt_s = cellfun(@(x) x.fit_rt_s, stat, 'uni', 0); fit_rt_s = cell2mat(fit_rt_s); fit_rt_mean_s = mean(fit_rt_s, 2); fit_rt_sem_s = std(fit_rt_s, 0, 2) / sqrt(nsubj);

% plot
fh = figure('color','w','Position',[100 100 350 150]);
subplot(1,2,1);
hold on;
h1 = plot(ucoh*100, fit_p_mean_ns, 'color', .5+.5*[0 0 0], 'linew', 1.5);
h2 = plot(ucoh*100, p_mean_ns, '.k', 'markers' ,7);
cerrorbar(ucoh*100, p_mean_ns, p_sem_ns, 'k');
plot(ucoh*100, fit_p_mean_s, 'color', .5+.5*[1 0 0], 'linew', 1.5);
plot(ucoh*100, p_mean_s, '.r', 'markers' ,7);
cerrorbar(ucoh*100, p_mean_s, p_sem_s, 'r');
axis square;
legend([h2,h1], {'Data', 'Model'}, 'location', 'southeast', 'box', 'off');
xlabel({'Stimulus strength','(% Morph)'})
ylabel(sprintf('P (response=%d)', opt.pos_cho))
ylim([0 1])
yticks(0:.25:1)

subplot(1,2,2);
hold on;
h1 = plot(ucoh*100, fit_rt_mean_ns, 'color', .5+.5*[0 0 0], 'linew', 1.5);
plot(ucoh*100, rtm_mean_ns, '.k', 'markers' ,7);
cerrorbar(ucoh*100, rtm_mean_ns, rtm_sem_ns, 'k');
h2 = plot(ucoh*100, fit_rt_mean_s, 'color', .5+.5*[1 0 0], 'linew', 1.5);
plot(ucoh*100, rtm_mean_s, '.r', 'markers' ,7);
cerrorbar(ucoh*100, rtm_mean_s, rtm_sem_s, 'r');
axis square;
legend([h1,h2], {'Non-switch', 'Switch'}, 'location', 'northwest', 'box', 'off');
xlabel({'Stimulus strength','(% Morph)'})
ylabel('Reaction time (s)')
ylim(opt.ylim)

end


function [fh, stat] = show_fit_choiceRT(coh_ns, resp_ns, fit_pchoice1_ns, rt_ns, fit_rt_ns, ...
    coh_s, resp_s, fit_pchoice1_s, rt_s, fit_rt_s, ...
    opt)

ucoh = unique(coh_ns);

% choice
[p_ns, pse_ns] = calcGroupMean(resp_ns==opt.pos_cho, coh_ns, ucoh, 'binary');
if opt.pos_cho==1
    fit_p_ns = calcGroupMean(fit_pchoice1_ns, coh_ns, ucoh);
else
    fit_p_ns = calcGroupMean(1 - fit_pchoice1_ns, coh_ns, ucoh);
end

[p_s, pse_s] = calcGroupMean(resp_s==opt.pos_cho, coh_s, ucoh, 'binary');
if opt.pos_cho==1
    fit_p_s = calcGroupMean(fit_pchoice1_s, coh_s, ucoh);
else
    fit_p_s = calcGroupMean(1 - fit_pchoice1_s, coh_s, ucoh);
end

% rt
[~, rtse_ns] = calcGroupMean(rt_ns, coh_ns, ucoh);
switch opt.RT_stat
    case 'mean'
        rtm_ns = calcGroupMean(rt_ns, coh_ns, ucoh);
        fit_rt_ns = calcGroupMean(fit_rt_ns, coh_ns, ucoh);
end

[~, rtse_s] = calcGroupMean(rt_s, coh_s, ucoh);
switch opt.RT_stat
    case 'mean'
        rtm_s = calcGroupMean(rt_s, coh_s, ucoh);
        fit_rt_s = calcGroupMean(fit_rt_s, coh_s, ucoh);
end

% plot
if opt.plot_idv_fh
    fh = figure('color', 'w', 'position', [100 100 250 150]);
    subplot(1,2,1);
    hold on;
    h1 = plot(ucoh*100, fit_p_ns, 'color', .5+.5*[0 0 0], 'linew', 1.5);
    h2 = plot(ucoh*100, p_ns, '.k', 'markers' ,7);
    cerrorbar(ucoh*100, p_ns, pse_ns, 'k');
    plot(ucoh*100, fit_p_s, 'color', .5+.5*[1 0 0], 'linew', 1.5);
    plot(ucoh*100, p_s, '.r', 'markers' ,7);
    cerrorbar(ucoh*100, p_s, pse_s, 'r');
    title(['subj=',opt.idv_subj]);
    axis square;
    % legend([h2,h1], {'Data', 'Model'}, 'location', 'southeast', 'box', 'off');
    xlabel({'Stimulus strength','(% Morph)'})
    ylabel(sprintf('P (response=%d)', opt.pos_cho))
    ylim([0 1])
    yticks(0:.25:1)
    
    subplot(1,2,2);
    hold on;
    plot(ucoh*100, fit_rt_ns, 'color', .5+.5*[0 0 0], 'linew', 1.5);
    plot(ucoh*100, rtm_ns, '.k', 'markers' ,7);
    cerrorbar(ucoh*100, rtm_ns, rtse_ns, 'k');
    plot(ucoh*100, fit_rt_s, 'color', .5+.5*[1 0 0], 'linew', 1.5);
    plot(ucoh*100, rtm_s, '.r', 'markers' ,7);
    cerrorbar(ucoh*100, rtm_s, rtse_s, 'r');
    axis square;
    xlabel({'Stimulus strength','(% Morph)'})
    ylabel('Reaction time (s)')
else
    fh = [];
end

% save
stat = struct('ucoh', ucoh, ...
    'p_ns', p_ns, 'fit_p_ns', fit_p_ns, 'p_s', p_s, 'fit_p_s', fit_p_s, ...
    'rtm_ns', rtm_ns, 'fit_rt_ns', fit_rt_ns, 'rtm_s', rtm_s, 'fit_rt_s', fit_rt_s);

end