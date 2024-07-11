function fh = run_show_fit_drt_csi(sim_file, opt)

opt.color = {[0 160 233]/255, [243 152 0]/255};
opt.sample = 'linear';

% plot
stat = load(sim_file).stat;
fh = show_drt_csi_average(stat, opt);

end



function fh = show_drt_csi_average(stat, opt)

% data averaging
nsubj = length(stat);
ucsi = stat{1}.ucsi;
ucsi_model = stat{1}.ucsi_model;
beta_0s = cellfun(@(x) x.beta_0s, stat, 'uni', 0); beta_0s = cell2mat(beta_0s); beta_0s_mean = mean(beta_0s, 2); beta_0s_sem = std(beta_0s, 0, 2) / sqrt(nsubj);
beta_0s_prep = cellfun(@(x) x.beta_0s_prep, stat, 'uni', 0); beta_0s_prep = cell2mat(beta_0s_prep); beta_0s_prep_mean = mean(beta_0s_prep, 2); beta_0s_prep_sem = std(beta_0s_prep, 0, 2) / sqrt(nsubj);
beta_0s_ktkr = cellfun(@(x) x.beta_0s_main, stat, 'uni', 0); beta_0s_ktkr = cell2mat(beta_0s_ktkr); beta_0s_ktkr_mean = mean(beta_0s_ktkr, 2); beta_0s_ktkr_sem = std(beta_0s_ktkr, 0, 2) / sqrt(nsubj);

% plot
fh = figure('color','w','Position',[100 100 350 150]);
subplot(1,2,1);
hold on
if strcmp(opt.sample,'linear')
    plot(ucsi, beta_0s_mean, '.', 'markers' ,7, 'Color', 'k'); % data point
    cerrorbar(ucsi, beta_0s_mean, beta_0s_sem, 'Color', 'k'); % error bar
    h1 = plot(ucsi_model, beta_0s_prep_mean, 'Color', opt.color{1}, 'LineWidth', 1.5);
    h2 = plot(ucsi_model, beta_0s_ktkr_mean, 'Color', opt.color{2}, 'LineWidth', 1.5);
    xlabel('CSI (s)')
    ylabel('\DeltaRT (s)')
    ylim(opt.ylim)
    plot(xlim, [0 0], '--', 'Color', [.5 .5 .5]);
    legend([h2, h1], {'Main model', 'Prep model'}, 'location', 'best', 'box', 'off');
end

end



