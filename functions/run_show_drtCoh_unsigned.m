function [fh, fh_idv, stat] = run_show_drtCoh_unsigned(cond_switch, coh, rt, subj, opt)

opt.log = true;
opt.plot_idv_fh = 0;

nsubj = length(opt.subj_list);
for s = 1:nsubj
    I = strcmp(subj, opt.subj_list{s});
    [stat{s}, fh_idv{s}] = show_choiceRT_unsigned(cond_switch(I), coh(I), rt(I), opt);
end
fh = show_choiceRT_unsigned_average(stat, opt);

end

function fh = show_choiceRT_unsigned_average(stat, opt)

nsubj = length(stat);
ucoh = stat{1}.ucoh;
lucoh = stat{1}.lucoh;

%% data averaging
% drt: calculate mean and sem across subjects
drt = cellfun(@(x) x.drt, stat, 'uni', 0); drt = cell2mat(drt); drt_mean = mean(drt, 2); drt_sem = std(drt, 0, 2) / sqrt(nsubj);

%% plot
fh = figure('color','w','Position',[100 100 350 150]);
subplot(1,2,1)
hold on
plot(lucoh, drt_mean, '.', 'markers' ,7, 'Color', 'k');
e = errorbar(lucoh, drt_mean, drt_sem, 'Color', 'k', 'LineStyle', 'none'); % errorbar
e.CapSize = 0;
axis square;
ucohl = arrayfun(@(x)num2str(x, '%g'), ucoh*100, 'uni', 0);
ucohl(2:2:end) = {''};
xticks(lucoh)
xticklabels(ucohl)
xlabel(opt.xlabel)
ylabel('\DeltaRT (s)')
ylim(opt.ylim)

end


function [stat, fh] = show_choiceRT_unsigned(cond, coh, rt, opt)

% coh, ucoh, lucoh and fcoh
coh = abs(coh);
ucoh = unique(coh);
if opt.log
    lucoh = log(ucoh);
    mincoh = ucoh(1);
    if mincoh == 0
        mincoh = ucoh(2)/ucoh(3) * ucoh(2);
        lucoh(1) = log(mincoh);
    end
else
    lucoh = ucoh;
    mincoh = 0;
end

% data averaging (rt)
rt1 = calcGroupMean(rt(cond==1), coh(cond==1), ucoh); % p1 ..averaged rt under cond1; pse1 ..standard error
rt2 = calcGroupMean(rt(cond==2), coh(cond==2), ucoh);
drt = rt2 - rt1;

% plot
if opt.plot_idv_fh
    fh = figure('color','w','Position',[100 100 200 100]);
    subplot(1,2,1)
    hold on
    plot(lucoh, drt, '.', 'markers' ,7, 'Color', 'k');
    axis square;
    ucohl = arrayfun(@(x)num2str(x, '%g'), ucoh*100, 'uni', 0);
    ucohl(2:2:end) = {''};
    xticks(lucoh)
    xticklabels(ucohl)
    xlabel(opt.xlabel)
    ylabel('\DeltaRT (s)')
else
    fh = [];
end

% save
stat = struct('ucoh', ucoh, 'lucoh', lucoh, 'drt', drt);

end

