function fh = run_show_fit_kernel(sim_file, opt)

opt.plot_idv_fh = false;

stat = load(sim_file).stat;
fh = show_kernel_average(stat, opt);

end


function fh = show_kernel_average(stat, opt)

nfeature = length(stat{1}.stkernel1);
t = stat{1}.t;
t_sim = stat{1}.t_sim;

stkernel1 = cellfun(@(x) x.stkernel1, stat, 'uni', 0); [stkernel1_mean, stkernel1_sem] = calKernelMean(stkernel1);
rekernel1 = cellfun(@(x) x.rekernel1, stat, 'uni', 0); [rekernel1_mean, rekernel1_sem] = calKernelMean(rekernel1);
stkernel2 = cellfun(@(x) x.stkernel2, stat, 'uni', 0); [stkernel2_mean, stkernel2_sem] = calKernelMean(stkernel2);
rekernel2 = cellfun(@(x) x.rekernel2, stat, 'uni', 0); [rekernel2_mean, rekernel2_sem] = calKernelMean(rekernel2);
stkernel_sim = cellfun(@(x) x.stkernel_sim, stat, 'uni', 0); [stkernel_sim_mean, stkernel_sim_sem] = calKernelMean(stkernel_sim);
rekernel_sim = cellfun(@(x) x.rekernel_sim, stat, 'uni', 0); [rekernel_sim_mean, rekernel_sim_sem] = calKernelMean(rekernel_sim);

fh = figure('color', 'w', 'position', [100 100 350 150*nfeature]);
ax = nan(nfeature,2);
for f = 1:nfeature
    ax(f,1) = subplot(nfeature,2,2*f-1);
    hold on;
    if opt.filter_cond_switch==1; plot_trace(t, stkernel1_mean{f}*100, stkernel1_sem{f}*100, [0 0 0]); end
    if opt.filter_cond_switch==2; plot_trace(t, stkernel2_mean{f}*100, stkernel2_sem{f}*100, [1 0 0]); end
    plot(t_sim, stkernel_sim_mean{f}*100, '-', 'color', [.3 .3 .3], 'linew', 3);
    axis square;
    plot(xlim, [0 0], '--', 'color', [.5 .5 .5]);
    xlabel('From stimulus on (s)');
    ylim(opt.ylim)
    ylabel('Kernel (% Morph)')
    
    ax(f,2) = subplot(nfeature,2,2*f);
    hold on;
    if opt.filter_cond_switch==1; plot_trace(fliplr(-t), rekernel1_mean{f}*100, rekernel1_sem{f}*100, [0 0 0]); end
    if opt.filter_cond_switch==2; plot_trace(fliplr(-t), rekernel2_mean{f}*100, rekernel2_sem{f}*100, [1 0 0]); end
    plot(fliplr(-t_sim), rekernel_sim_mean{f}*100, '-', 'color', [.3 .3 .3], 'linew', 3);
    axis square;
    plot(xlim, [0 0], '--', 'color', [.5 .5 .5]);
    xlabel('From response (s)');
    set(gca, 'YAxisLocation', 'right');
    ylim(opt.ylim)
end

end


function lh = plot_trace(x, y, se, color)

hold on;

if isempty(se)
    se = zeros(size(y));
end

palecol = color * .5 + [1 1 1] * .5;
x2 = [x(:); flipud(x(:))];
y2 = [y(:)-se(:); flipud(y(:)+se(:))];
h = fill(x2,y2,palecol);
set(h, 'EdgeColor', palecol, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% draw main trace
lh = plot(x,y,'Color', color, 'LineWidth', 0.5);

end


function [kernel_mean, kernel_sem] = calKernelMean(kernel)

nsubj = length(kernel);
nfeature = length(kernel{1});
medRT = min(cellfun(@(x) length(x{1}), kernel));

kernel_mean = cell(nfeature,1);
kernel_sem  = cell(nfeature,1);
for f = 1:nfeature
    kernel_allSubj = nan(nsubj, medRT);
    for s = 1:nsubj
        kernel_allSubj(s,:) = kernel{s}{f};
    end
    kernel_mean{f} = mean(kernel_allSubj,1);
    kernel_sem{f}  = std (kernel_allSubj,0,1)/sqrt(nsubj);
end

end

