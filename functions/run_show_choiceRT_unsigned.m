function [fh, fh_idv, stat] = run_show_choiceRT_unsigned(cond_switch, coh, resp, targ_cor, rt, subj, opt)

opt.plot_idv_fh = false;
opt.log = true;
opt.constant_off = true;


nsubj = length(opt.subj_list);
for s = 1:nsubj
    I = strcmp(subj, opt.subj_list{s});
    opt.subj_idx = s;
    [stat{s}, fh_idv{s}] = show_choiceRT_unsigned(cond_switch(I), coh(I), resp(I), targ_cor(I), rt(I), opt);
end
fh = show_choiceRT_unsigned_average(stat, opt);

end

function fh = show_choiceRT_unsigned_average(stat, opt)

nsubj = length(stat);
ucoh = stat{1}.ucoh;
lucoh = stat{1}.lucoh;
fcoh = stat{1}.fcoh;

%% data averaging
% resp: calculate mean and sem across subjects
resp1 = cellfun(@(x) x.resp1, stat, 'uni', 0); resp1 = cell2mat(resp1); p1_mean = mean(resp1, 2); p1_sem = std(resp1, 0, 2) / sqrt(nsubj);
resp2 = cellfun(@(x) x.resp2, stat, 'uni', 0); resp2 = cell2mat(resp2); p2_mean = mean(resp2, 2); p2_sem = std(resp2, 0, 2) / sqrt(nsubj);
% rt: calculate mean and sem across subjects
rt1 = cellfun(@(x) x.rt1, stat, 'uni', 0); rt1 = cell2mat(rt1); rt1_mean = mean(rt1, 2); rt1_sem = std(rt1, 0, 2) / sqrt(nsubj);
rt2 = cellfun(@(x) x.rt2, stat, 'uni', 0); rt2 = cell2mat(rt2); rt2_mean = mean(rt2, 2); rt2_sem = std(rt2, 0, 2) / sqrt(nsubj);

%% fit curve averaging
% resp
fresp1 = cellfun(@(x) x.fresp1, stat, 'uni', 0); fresp1 = cell2mat(fresp1); fresp1 = mean(fresp1, 2);
fresp2 = cellfun(@(x) x.fresp2, stat, 'uni', 0); fresp2 = cell2mat(fresp2); fresp2 = mean(fresp2, 2);
% rt
frt1 = cellfun(@(x) x.frt1, stat, 'uni', 0); frt1 = cell2mat(frt1); frt1 = mean(frt1, 2);
frt2 = cellfun(@(x) x.frt2, stat, 'uni', 0); frt2 = cell2mat(frt2); frt2 = mean(frt2, 2);

%% plot
fh = figure('color','w','Position',[100 100 350 150]);
subplot(1,2,1)
hold on
if opt.log
    plot(log(fcoh), fresp1, 'Color', 'k'); % fit curve
    plot(log(fcoh), fresp2, 'Color', 'r');
else
    plot(fcoh, fresp1, 'Color', 'k');
    plot(fcoh, fresp2, 'Color', 'r');
end
plot(lucoh, p1_mean, '.', 'markers' ,7, 'Color', 'k');
plot(lucoh, p2_mean, '.', 'markers' ,7, 'Color', 'r');
cerrorbar(lucoh, p1_mean, p1_sem, 'Color', 'k');
cerrorbar(lucoh, p2_mean, p2_sem, 'Color', 'r');
axis square;
ucohl = arrayfun(@(x)num2str(x, '%g'), ucoh*100, 'uni', 0);
ucohl(2:2:end) = {''};
ylim([.5 1])
yticks(linspace(.5,1,5))
yticklabels({'0.5','','0.75','','1'})
xticks(lucoh)
xticklabels(ucohl)
xlabel({'Stimulus strength','(% Morph)'})
ylabel('P (correct)')

subplot(1,2,2)
hold on;
plot(log(fcoh), frt1, 'Color', 'k');
plot(lucoh, rt1_mean, '.', 'markers' ,7, 'Color', 'k');
cerrorbar(lucoh, rt1_mean, rt1_sem, 'Color', 'k');
plot(log(fcoh), frt2, 'Color', 'r');
plot(lucoh, rt2_mean, '.', 'markers' ,7, 'Color', 'r');
cerrorbar(lucoh, rt2_mean, rt2_sem, 'Color', 'r');
axis square;
ucohl = arrayfun(@(x)num2str(x, '%g'), ucoh*100, 'uni', 0);
ucohl(2:2:end) = {''};
yticks(linspace(.2,.4,5))
yticklabels({'0.2','','0.3','','0.4'})
xticks(lucoh)
xticklabels(ucohl)
xlabel({'Stimulus strength','(% Morph)'})
ylabel('Reaction time (s)')

end


function [stat, fh] = show_choiceRT_unsigned(cond, coh, resp, targ_cor, rt, opt)

% corr
corr = resp==targ_cor;

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
fcoh = linspace(mincoh, max(coh), 100)';

% data averaging (choice)
[resp1, respse1] = calcGroupMean(corr(cond==1), coh(cond==1), ucoh, 'binary'); % p1 ..prob under cond1; pse1 ..standard error
[resp2, respse2] = calcGroupMean(corr(cond==2), coh(cond==2), ucoh, 'binary');

% psychometric function
if opt.constant_off
    [alpha, ~, glmfit_stat_alpha] = glmfit([coh, coh.*(cond==2)], corr, 'binomial', 'link', 'logit', 'Constant', 'off'); % logit(P) = (a1 + a1s * cond) * coh
    fresp1 = glmval(alpha, [fcoh, fcoh.*zeros(size(fcoh))], 'logit', 'Constant', 'off'); % fresp1 ..function_resp_1
    fresp2 = glmval(alpha, [fcoh, fcoh.*ones(size(fcoh))], 'logit', 'Constant', 'off');
else
    [alpha, ~, glmfit_stat_alpha] = glmfit([cond==2, coh, coh.*(cond==2)], corr, 'binomial', 'link', 'logit'); % logit(P) = a0 + a0s * cond + (a1 + a1s * cond) * coh
    fresp1 = glmval(alpha, [zeros(size(fcoh)), fcoh, fcoh.*zeros(size(fcoh))], 'logit'); % fresp1 ..function_resp_1
    fresp2 = glmval(alpha, [ones(size(fcoh)), fcoh, fcoh.*ones(size(fcoh))], 'logit');
end

% data averaging (rt)
[rt1, rtse1] = calcGroupMean(rt(cond==1), coh(cond==1), ucoh); % p1 ..averaged rt under cond1; pse1 ..standard error
[rt2, rtse2] = calcGroupMean(rt(cond==2), coh(cond==2), ucoh);

% chronometric function
beta = nlinfit([cond==2, coh], rt, @hyperbolic_fun, 0.1*ones(1,4));
frt1 = hyperbolic_fun(beta, [zeros(size(fcoh)), abs(fcoh)]);
frt2 = hyperbolic_fun(beta, [ones(size(fcoh)), abs(fcoh)]);

% plot
if opt.plot_idv_fh
    fh = figure('color','w','Position',[200+200*opt.subj_idx 100 200 100]);
    subplot(1,2,1)
    hold on
    if opt.log
        plot(log(fcoh), fresp1, 'Color', 'k'); % fit curve
        plot(log(fcoh), fresp2, 'Color', 'r');
    else
        plot(fcoh, fresp1, 'Color', 'k');
        plot(fcoh, fresp2, 'Color', 'r');
    end
    plot(lucoh, resp1, '.', 'markers' ,7, 'Color', 'k');
    plot(lucoh, resp2, '.', 'markers' ,7, 'Color', 'r');
    cerrorbar(lucoh, resp1, respse1, 'Color', 'k');
    cerrorbar(lucoh, resp2, respse2, 'Color', 'r');
    axis square;
    ucohl = arrayfun(@(x)num2str(x, '%g'), ucoh*100, 'uni', 0);
    ucohl(2:2:end) = {''};
    title(['subj ', opt.subj_list{opt.subj_idx}]);
    
    subplot(1,2,2)
    hold on;
    plot(log(fcoh), frt1, 'Color', 'k');
    plot(lucoh, rt1, '.', 'markers' ,7, 'Color', 'k');
    cerrorbar(lucoh, rt1, rtse1, 'Color', 'k');
    plot(log(fcoh), frt2, 'Color', 'r');
    plot(lucoh, rt2, '.', 'markers' ,7, 'Color', 'r');
    cerrorbar(lucoh, rt2, rtse2, 'Color', 'r');
    axis square;
    ucohl = arrayfun(@(x)num2str(x, '%g'), ucoh*100, 'uni', 0);
    ucohl(2:2:end) = {''};
else
    fh = [];
end

% save
stat = struct('ucoh', ucoh, 'lucoh', lucoh, 'resp1', resp1, 'resp2', resp2, 'rt1', rt1, 'rt2', rt2, ...
    'fcoh', fcoh, 'fresp1', fresp1, 'fresp2', fresp2, 'frt1', frt1, 'frt2', frt2, ...
    'alpha', alpha(:), 'beta', beta(:), ...
    'glmfit_stat_alpha', glmfit_stat_alpha, ...
    'constant_off', opt.constant_off);

end


function y = hyperbolic_fun(beta, X)

% X(:,1) ..cond
% X(:,2) ..coh

cond = X(:,1);
coh = X(:,2);

y = beta(1) + beta(2) .* cond + ...
    beta(3) ./ coh .* tanh(beta(4) .* coh);
I = coh==0; % coh=0
y(I) = beta(1) + beta(2) .* cond(I) + ...
    beta(3) * beta(4);

end