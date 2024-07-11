function [fh, fh_idv, stat] = run_show_choiceRT_mot_col(cond_switch, coh, resp, rt, subj, opt)

opt.plot_idv_fh = false;

nsubj = length(opt.subj_list);
for s = 1:nsubj
    I = strcmp(subj, opt.subj_list{s});
    opt.subj_idx = s;
    [stat{s}, fh_idv{s}] = show_choiceRT(cond_switch(I), coh(I), resp(I), rt(I), opt);
end
fh = show_choiceRT_average(stat, opt);

end

function fh = show_choiceRT_average(stat, opt)

nsubj = length(stat);
ucoh = stat{1}.ucoh;
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
subplot(1,2,1);
hold on;
plot(fcoh, fresp1, 'Color', opt.color{1}); % averaged regression curve under cond 1
plot(ucoh, p1_mean, '.', 'markers' ,7, 'Color', opt.color{1}); % data points under cond 1
cerrorbar(ucoh, p1_mean, p1_sem, 'Color', opt.color{1}); % errorbar
plot(fcoh, fresp2, 'Color', opt.color{2}); % averaged regression curve under cond 2
plot(ucoh, p2_mean, '.', 'markers' ,7, 'Color', opt.color{2}); % data points under cond 2
cerrorbar(ucoh, p2_mean, p2_sem, 'Color', opt.color{2}); % errorbar
axis square;
ylim([0 1]); yticks(0:.25:1); ylabel('P (response=2)')
xlim(opt.xlim); xlabel(opt.xlabel)

subplot(1,2,2);
hold on;
h1 = plot(fcoh, frt1, 'Color', opt.color{1}); % averaged regression curve under cond 1
plot(ucoh, rt1_mean, '.', 'markers' ,7, 'Color', opt.color{1}); % data points under cond 1
cerrorbar(ucoh, rt1_mean, rt1_sem, 'Color', opt.color{1}); % errorbar
h2 = plot(fcoh, frt2, 'Color', opt.color{2}); % averaged regression curve under cond 2
plot(ucoh, rt2_mean, '.', 'markers' ,7, 'Color', opt.color{2}); % data points under cond 2
cerrorbar(ucoh, rt2_mean, rt2_sem, 'Color', opt.color{2}); % errorbar
axis square;
legend([h1, h2], opt.legend, 'Location', 'best'); legend boxoff;
ylim(opt.ylim); ylabel('Reaction time (s)')
xlim(opt.xlim); xlabel(opt.xlabel)

end

function [stat, fh] = show_choiceRT(cond, coh, resp, rt, opt)

% data averaging (choice)
ucoh = unique(coh); % ucoh ..unique coh
[resp1, respse1] = calcGroupMean(resp(cond==1)==2, coh(cond==1), ucoh, 'binary'); % p1 ..prob under cond1; pse1 ..standard error
[resp2, respse2] = calcGroupMean(resp(cond==2)==2, coh(cond==2), ucoh, 'binary');

% psychometric function
fcoh = linspace(-1, 1, 100)'; % fcoh ..function_coh
alpha = glmfit([cond==2, coh, coh.*(cond==2)], resp==2, 'binomial', 'link', 'logit'); % logit(P) = a0 + a0s * cond + (a1 + a1s * cond) * coh
fresp1 = glmval(alpha, [zeros(size(fcoh)), fcoh, fcoh.*zeros(size(fcoh))], 'logit'); % fresp1 ..function_resp_1
fresp2 = glmval(alpha, [ones(size(fcoh)),  fcoh, fcoh.*ones(size(fcoh))], 'logit');

% data averaging (rt)
[rt1, rtse1] = calcGroupMean(rt(cond==1), coh(cond==1), ucoh); % p1 ..averaged rt under cond1; pse1 ..standard error
[rt2, rtse2] = calcGroupMean(rt(cond==2), coh(cond==2), ucoh);

% chronometric function
beta = nlinfit([cond==2, abs(coh)], rt, @hyperbolic_fun, 0.1*ones(1,4));
frt1 = hyperbolic_fun(beta, [zeros(size(fcoh)), abs(fcoh)]);
frt2 = hyperbolic_fun(beta, [ones(size(fcoh)), abs(fcoh)]);

% plot
if opt.plot_idv_fh
    fh = figure('color','w','Position',[200+200*opt.subj_idx 100 200 100]);
    subplot(1,2,1)
    hold on
    h1 = plot(fcoh, fresp1, 'Color', opt.color{1}); % regression curve
    plot(ucoh, resp1, '.', 'markers' ,7, 'Color', opt.color{1}); % data points
    cerrorbar(ucoh, resp1, respse1, 'Color', opt.color{1}); % error bar
    h2 = plot(fcoh, fresp2, 'Color', opt.color{2});
    plot(ucoh, resp2, '.', 'markers' ,7, 'Color', opt.color{2});
    cerrorbar(ucoh, resp2, respse2, 'Color', opt.color{2});
    axis square;
    legend([h1, h2], {'Nonswitch', 'Switch'}, 'Location', 'best'); legend boxoff; legend off
    title(['subj ', opt.subj_list{opt.subj_idx}]);
    ylim([0 1])
    yticks(0:.25:1)
    xlabel({'Stimulus strength','(% Morph)'})
    ylabel('P (response=2)')
    
    subplot(1,2,2)
    hold on
    plot(fcoh, frt1, 'Color', opt.color{1}); % fit curve
    plot(ucoh, rt1, '.', 'markers' ,7, 'Color', opt.color{1}); % data points
    cerrorbar(ucoh, rt1, rtse1, 'Color', opt.color{1}); % errorbar
    plot(fcoh, frt2, 'Color', opt.color{2});
    plot(ucoh, rt2, '.', 'markers' ,7, 'Color', opt.color{2});
    cerrorbar(ucoh, rt2, rtse2, 'Color', opt.color{2});
    axis square;
    xlabel({'Stimulus strength','(% Morph)'})
    ylabel('Reaction time (s)')
else
    fh = [];
end
   
% save
stat = struct('ucoh', ucoh, 'resp1', resp1, 'resp2', resp2, 'rt1', rt1, 'rt2', rt2, ...
    'fcoh', fcoh, 'fresp1', fresp1, 'fresp2', fresp2, 'frt1', frt1, 'frt2', frt2, ...
    'alpha', alpha(:), 'beta', beta(:));

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
    beta(3) .* beta(4);
end