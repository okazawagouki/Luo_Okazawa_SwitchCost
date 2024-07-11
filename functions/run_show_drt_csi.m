function [fh, stat] = run_show_drt_csi(cond_switch, coh, csi, rt, subj, opt)

opt.plot_idv_fh = false;
opt.sample = 'linear';

% csi sampling
[opt.csi_edge, opt.ucsi] = sampling(csi, opt);
opt.fcsi = linspace(min(csi), max(csi), 100)';

% plot
nsubj = length(opt.subj_list);
for s = 1:nsubj
    I = strcmp(subj, opt.subj_list{s});
    stat{s} = show_drt_csi(cond_switch(I), coh(I), csi(I), rt(I), opt);
end
fh = show_drt_csi_average(stat, opt);

end

function fh = show_drt_csi_average(stat, opt)

nsubj = length(stat);

% data averaging
ucsi = stat{1}.ucsi;
beta_0s = cellfun(@(x) x.beta_0s, stat, 'uni', 0); beta_0s = cell2mat(beta_0s);
beta_0s_mean = mean(beta_0s, 2); beta_0s_sem = std(beta_0s, 0, 2) / sqrt(nsubj);

% fitting
fcsi = stat{1}.fcsi;
f_beta_0s = cellfun(@(x) x.f_beta_0s, stat, 'uni', 0); f_beta_0s = cell2mat(f_beta_0s);
f_beta_0s_mean = mean(f_beta_0s, 2); f_beta_0s_sem = std(f_beta_0s, 0, 2) / sqrt(nsubj);

% plot
fh = figure('color','w','Position',[100 100 350 150]);
subplot(1,2,1);
hold on
if strcmp(opt.sample,'log')
    % x-axis is log scale
    plot_trace(log(fcsi), f_beta_0s_mean, f_beta_0s_sem, [.5 .5 .5]);
    plot(log(ucsi), beta_0s_mean, '.', 'markers' ,7, 'Color', 'k');
    cerrorbar(log(ucsi), beta_0s_mean, beta_0s_sem, 'Color', 'k');
    xlabel('CSI (s)')
    ylabel('\DeltaRT (s)')
    ylim(opt.ylim)
    xlim(log([min(fcsi), max(fcsi)]))
    xticks(log(linspace(min(fcsi), max(fcsi), 5)))
    xticklabels(num2str(linspace(min(fcsi), max(fcsi), 5)','%.1f'))
elseif strcmp(opt.sample,'linear')
    % plot(fcsi, fdrt_mean, 'Color', 'k'); % only fit curve
    plot_trace(fcsi, f_beta_0s_mean, f_beta_0s_sem, [.5 .5 .5]); % fit curve with shade
    plot(ucsi, beta_0s_mean, '.', 'markers' ,7, 'Color', 'k'); % data point
    cerrorbar(ucsi, beta_0s_mean, beta_0s_sem, 'Color', 'k');
    xlabel('CSI (s)')
    ylabel('\DeltaRT (s)')
    ylim(opt.ylim)
    yticks(0:0.1:0.4);
    yticklabels({'0','','0.2','','0.4'})
end
ax = get(gcf,'Children');
axis normal
set(ax(end), 'Position', [0.2 0.3 0.45 0.5])
set(ax(end), 'TickLength', [0.015 0.015])

end

function [stat, fh] = show_drt_csi(cond_switch, coh, csi, rt, opt)

% fit each csi bin to chronometric function
for t = 1:length(opt.ucsi)
    if t==1
        I = csi>=opt.csi_edge(t) & csi<=opt.csi_edge(t+1);
    else
        I = csi> opt.csi_edge(t) & csi<=opt.csi_edge(t+1);
    end
    beta{t} = nlinfit([cond_switch(I)==2, abs(coh(I)), csi(I)], rt(I), @hyperbolic_fun, 0.1*ones(1,4));
end
beta_0s = cellfun(@(x) x(2), beta);

% linear regression of beta_0s and csi
b = regress(beta_0s(:), [ones(size(opt.ucsi(:))), opt.ucsi]);
f_beta_0s = b(1) + b(2)*opt.fcsi;

% plot
if opt.plot_idv_fh
    fh = figure('color','w','Position',[100 100 200 100]);
    subplot(1,2,1);
    hold on
    plot(opt.ucsi, beta_0s, '.', 'markers' ,7, 'Color', 'k'); % data points
    plot(opt.fcsi, f_beta_0s, 'Color', 'k'); % fit curve
else
    fh = [];
end

% save
% beta ..fit to chronometric function
% b ..fit to linear regression
stat = struct('ucsi', opt.ucsi, 'beta_0s', beta_0s(:), ...
    'fcsi', opt.fcsi, 'f_beta_0s', f_beta_0s, 'b', b);
stat.beta = beta(:);

end

function y = hyperbolic_fun(beta, X)
% X(:,1) ..cond
% X(:,2) ..coh
% X(:,3) ..t

cond = X(:,1);
coh = X(:,2);

y = beta(1) + beta(2) .* cond + ...
    beta(3) ./ coh .* tanh(beta(4) .* coh);
I = coh==0; % coh=0
y(I) = beta(1) + beta(2) .* cond(I) + ...
    beta(3) .* beta(4);
end