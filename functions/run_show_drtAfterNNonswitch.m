function [fh, fh_idv, stat] = run_show_drtAfterNNonswitch(cond_switch, n_after_nonswitch, rt, subj, opt)

opt.plot_idv_fh = false;

nsubj = length(opt.subj_list);
for s = 1:nsubj
    I = strcmp(subj, opt.subj_list{s});
    [stat{s}, fh_idv{s}] = show_drtAfterNSwitch(cond_switch(I), n_after_nonswitch(I), rt(I), opt);
end
fh = show_drtAfterNSwitch_average(stat, opt);

end

function fh = show_drtAfterNSwitch_average(stat, opt)

nsubj = length(stat);
u_n_after_nonswitch = stat{1}.u_n_after_nonswitch;

%% data averaging
% drt: calculate mean and sem across subjects
drt = cellfun(@(x) x.drt, stat, 'uni', 0); drt = cell2mat(drt); drt_mean = mean(drt, 2); drt_sem = std(drt, 0, 2) / sqrt(nsubj);
fdrt = cellfun(@(x) x.fdrt, stat, 'uni', 0); fdrt = cell2mat(fdrt); fdrt_mean = mean(fdrt, 2); fdrt_sem = std(fdrt, 0, 2) / sqrt(nsubj);

%% plot
fh = figure('color','w','Position',[100 100 350 150]);
subplot(1,2,1);
hold on;
% plot_trace(u_n_after_nonswitch, fdrt_mean, fdrt_sem, [.5 .5 .5]); % fit curve with shade
plot(u_n_after_nonswitch, drt_mean, '.', 'markers' ,7, 'Color', 'k');
cerrorbar(u_n_after_nonswitch, drt_mean, drt_sem, 'Color', 'k');
axis square;
xlabel({'Number of non-switch trials', 'before task switch'})
ylabel('\Delta RT (s)')
xlim(opt.xlim)
ylim(opt.ylim)
xticks([1 2 3 4])
plot(xlim, [0 0], '--', 'Color', [.5 .5 .5])

end

function [stat, fh] = show_drtAfterNSwitch(cond_switch, n_after_nonswitch, rt, opt)

% data averaging
u_n_after_nonswitch = unique(n_after_nonswitch(~isnan(n_after_nonswitch))); % 1 ..after 1 non-switch trial
[rt_s] = calcGroupMean(rt(cond_switch==2), n_after_nonswitch(cond_switch==2), u_n_after_nonswitch);
meanrt_ns = mean(rt(cond_switch==1));
drt = rt_s - meanrt_ns;

% cut off and regression
u_n_after_nonswitch = u_n_after_nonswitch(1:4);
drt = drt(1:4);
b = regress(drt, [ones(length(u_n_after_nonswitch),1), u_n_after_nonswitch]);
fdrt = [ones(length(u_n_after_nonswitch),1), u_n_after_nonswitch] * b;

% plot
if opt.plot_idv_fh
    fh = figure('color','w','Position',[100 100 200 100]);
    subplot(1,2,1)
    hold on
    plot(u_n_after_nonswitch, drt, '.', 'markers' ,7, 'Color', 'k');
    plot(u_n_after_nonswitch, fdrt, 'k');
    axis square;
    xlabel({'Number of non-switch trials', 'before task switch'})
    ylabel('\Delta RT (s)')
    ylim(opt.ylim)
else
    fh = [];
end

% save
stat = struct('u_n_after_nonswitch', u_n_after_nonswitch, 'drt', drt, ...
    'b', b, 'fdrt', fdrt);

end
