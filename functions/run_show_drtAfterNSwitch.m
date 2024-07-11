function [fh, fh_idv, stat] = run_show_drtAfterNSwitch(cond_switch, n_after_switch, rt, subj, opt)

opt.plot_idv_fh = false;

nsubj = length(opt.subj_list);
for s = 1:nsubj
    I = strcmp(subj, opt.subj_list{s});
    [stat{s}, fh_idv{s}] = show_drtAfterNSwitch(cond_switch(I), n_after_switch(I), rt(I), opt);
end
fh = show_drtAfterNSwitch_average(stat, opt);

end

function fh = show_drtAfterNSwitch_average(stat, opt)

nsubj = length(stat);
u_n_after_switch = stat{1}.u_n_after_switch;

%% data averaging
% drt: calculate mean and sem across subjects
drt = cellfun(@(x) x.drt, stat, 'uni', 0); drt = cell2mat(drt); drt_mean = mean(drt, 2); drt_sem = std(drt, 0, 2) / sqrt(nsubj);

%% plot
fh = figure('color','w','Position',[100 100 350 150]);
subplot(1,2,1);
hold on;
h1 = plot(u_n_after_switch(2:end)+1, drt_mean(2:end), '.', 'markers' ,7, 'Color', 'k');
cerrorbar(u_n_after_switch(2:end)+1, drt_mean(2:end), drt_sem(2:end), 'Color', 'k');
h2 = plot(u_n_after_switch(1)+1, drt_mean(1), '.', 'markers' ,7, 'Color', 'r');
cerrorbar(u_n_after_switch(1)+1, drt_mean(1), drt_sem(1), 'Color', 'r');
axis square;
xlim(opt.xlim)
ylim(opt.ylim)
xlabel({'Number of trials', 'after task switch'})
ylabel('\DeltaRT (s)')
plot(xlim, [0 0], '--', 'Color', [.5 .5 .5])

end

function [stat, fh] = show_drtAfterNSwitch(cond_switch, n_after_switch, rt, opt)

% data averaging
u_n_after_switch = unique(n_after_switch); % 0 ..switch; 1 ..after 1 switch trial
[rt_ns] = calcGroupMean(rt(cond_switch==1), n_after_switch(cond_switch==1), u_n_after_switch); % p1 ..averaged rt under cond1; pse1 ..standard error
[rt_s] = calcGroupMean(rt(cond_switch==2), n_after_switch(cond_switch==2), u_n_after_switch);
% in the previous version of code, I took the mean of the grouped mean
meanrt_ns = nanmean(rt_ns); 

% in the current version, I directly calculated mean rt from all ns trials, this number should be more accurate
% meanrt_ns = mean(rt(cond_switch==1)); 

drt = [rt_s(1); rt_ns(2:end)] - meanrt_ns;

% plot
if opt.plot_idv_fh
    fh = figure('color','w','Position',[100 100 200 100]);
    subplot(1,2,1)
    hold on
    plot(u_n_after_switch+1, drt, '.', 'markers' ,7, 'Color', 'k');
    axis square;
    xlabel({'Number of trials', 'after task switch'})
    ylabel('\DeltaRT (s)')
else
    fh = [];
end

% save
stat = struct('u_n_after_switch', u_n_after_switch, 'drt', drt);

end

