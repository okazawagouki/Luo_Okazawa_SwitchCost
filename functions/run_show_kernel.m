
function run_show_kernel(cond_switch, morph_level, cond, fluc, resp, subj, opt)

opt.medRT = round(prctile(cellfun(@(x) size(x,1), fluc), 60));

[~, ~, coh_main] = get_coh(morph_level, cond);

for s = 1:length(opt.subj_list)
    param_ns = load(opt.fit_file{s}).fitresult.model_param.final;
    k = reshape(param_ns(1:6), [6,1]);
    I_subj = strcmp(subj, opt.subj_list{s});
    fluc(cond==1 & I_subj) = cellfun(@(x) x(:,1:3)*k(1:3)/sum(k(1:3)), fluc(cond==1 & I_subj), 'uni', 0);
    fluc(cond==2 & I_subj) = cellfun(@(x) x(:,4:6)*k(4:6)/sum(k(4:6)), fluc(cond==2 & I_subj), 'uni', 0);
end
coh = coh_main;

nsubj = length(opt.subj_list);
stat = cell(1, nsubj);
for s = 1:nsubj
    I = strcmp(subj, opt.subj_list{s});
    [~, stat{s}] = show_kernel(cond_switch(I), coh(I), fluc(I), resp(I), opt);
end

show_kernel_average(stat, opt);

end


function fh = show_kernel_average(stat, opt)

nfeature = length(stat{1}.stkernel1);
t = stat{1}.t;

stkernel1 = cellfun(@(x) x.stkernel1, stat, 'uni', 0); [stkernel1_mean, stkernel1_sem] = calKernelMean(stkernel1);
rekernel1 = cellfun(@(x) x.rekernel1, stat, 'uni', 0); [rekernel1_mean, rekernel1_sem] = calKernelMean(rekernel1);
stkernel2 = cellfun(@(x) x.stkernel2, stat, 'uni', 0); [stkernel2_mean, stkernel2_sem] = calKernelMean(stkernel2);
rekernel2 = cellfun(@(x) x.rekernel2, stat, 'uni', 0); [rekernel2_mean, rekernel2_sem] = calKernelMean(rekernel2);

fh = figure('color', 'w', 'position', [100 100 350 150*nfeature]);
ax = nan(nfeature,2);
for f = 1:nfeature
    ax(f,1) = subplot(nfeature,2,2*f-1);
    hold on;
    plot_trace(t, stkernel1_mean{f}*100, stkernel1_sem{f}*100, [0 0 0]);
    plot_trace(t, stkernel2_mean{f}*100, stkernel2_sem{f}*100, [1 0 0]);
    axis square;
    plot(xlim, [0 0], '--', 'color', [.5 .5 .5]);
    xlabel('From stimulus on (s)');
    ylabel('Kernel (% Morph)');
    
    ax(f,2) = subplot(nfeature,2,2*f);
    hold on;
    plot_trace(fliplr(-t), rekernel1_mean{f}*100, rekernel1_sem{f}*100, [0 0 0]);
    plot_trace(fliplr(-t), rekernel2_mean{f}*100, rekernel2_sem{f}*100, [1 0 0]);
    axis square;
    plot(xlim, [0 0], '--', 'color', [.5 .5 .5]);
    xlabel('From response (s)');
    set(gca, 'YAxisLocation', 'right');
end
ax = get(gcf,'Children');
for n = 1:length(ax)
    set(ax(n), 'ylim', opt.ylim);
end

end


function [fh, stat] = show_kernel(cond_switch, coh, fluc, resp, opt)

opt.max_coh = 0.12; % maximum coherence of trials used for kernel analysis
opt.smoothing = 3; % smoothing of kernel, e.g. fspecial('average', [1, 3])
opt.fluc_t = 8/75; % duration of each fluctuation (sec)
opt.positive_cat = 2; % category (choice) of positive coherence
opt.plot_idv_fh = 0;

if ~isempty(opt.max_coh)
    ind = abs(coh) <= opt.max_coh;
    cond_switch = cond_switch(ind);
    coh = coh(ind);
    fluc = fluc(ind);
    resp = resp(ind);
end

ntrial = length(fluc);
medRT = opt.medRT;
nfeature = size(fluc{1},2);

stfluc = repmat({nan(ntrial, medRT)}, [nfeature,1]);
refluc = repmat({nan(ntrial, medRT)}, [nfeature,1]);
for n = 1:ntrial
    ln = min(size(fluc{n},1), medRT);
    if opt.smoothing>1
        fluc{n} = nanconv(fluc{n}, fspecial('average', [opt.smoothing, 1]), 'same');
    end
    for f = 1:nfeature
        stfluc{f}(n, 1:ln)         = fluc{n}(1:ln        , f) - coh(n);
        refluc{f}(n, end-ln+1:end) = fluc{n}(end-ln+1:end, f) - coh(n);
    end
end

stkernel1 = cell(nfeature,1); stkernel1_se = cell(nfeature,1); stkernel2 = cell(nfeature,1); stkernel2_se = cell(nfeature,1);
rekernel1 = cell(nfeature,1); rekernel1_se = cell(nfeature,1); rekernel2 = cell(nfeature,1); rekernel2_se = cell(nfeature,1);
for f = 1:nfeature
    I_ns = cond_switch==1; I_s = cond_switch==2;
    % non-switch
    stkernel1{f} = nanmean(stfluc{f}(resp==opt.positive_cat & I_ns,:)) - nanmean(stfluc{f}(resp==3-opt.positive_cat & I_ns,:));
    stkernel1_se{f} = calc_diff_SE(stfluc{f}(resp==1 & I_ns,:), stfluc{f}(resp==2 & I_ns,:));
    rekernel1{f} = nanmean(refluc{f}(resp==opt.positive_cat & I_ns,:)) - nanmean(refluc{f}(resp==3-opt.positive_cat & I_ns,:));
    rekernel1_se{f} = calc_diff_SE(refluc{f}(resp==1 & I_ns,:), refluc{f}(resp==2 & I_ns,:));
    % switch
    stkernel2{f} = nanmean(stfluc{f}(resp==opt.positive_cat & I_s,:)) - nanmean(stfluc{f}(resp==3-opt.positive_cat & I_s,:));
    stkernel2_se{f} = calc_diff_SE(stfluc{f}(resp==1 & I_s,:), stfluc{f}(resp==2 & I_s,:));
    rekernel2{f} = nanmean(refluc{f}(resp==opt.positive_cat & I_s,:)) - nanmean(refluc{f}(resp==3-opt.positive_cat & I_s,:));
    rekernel2_se{f} = calc_diff_SE(refluc{f}(resp==1 & I_s,:), refluc{f}(resp==2 & I_s,:));
end
t = (0:medRT-1) * opt.fluc_t;

if opt.plot_idv_fh
    fh = figure('color', 'w', 'position', [100 100 350 150*nfeature]);
    ax = nan(nfeature,2);
    for f = 1:nfeature
        ax(f,1) = subplot(nfeature,2,2*f-1);
        hold on;
        plot_trace(t, stkernel1{f}*100, stkernel1_se{f}*100, [0 0 0]);
        plot_trace(t, stkernel2{f}*100, stkernel2_se{f}*100, [1 0 0]);
        axis square;
        plot(xlim, [0 0], '--', 'color', [.5 .5 .5]);
        xlabel('From stimulus on (s)');
        ylabel('Kernel (% Morph)');
        ylim1 = get(gca,'ylim');
        
        ax(f,2) = subplot(nfeature,2,2*f);
        hold on;
        plot_trace(fliplr(-t), rekernel1{f}*100, rekernel1_se{f}*100, [0 0 0]);
        plot_trace(fliplr(-t), rekernel2{f}*100, rekernel2_se{f}*100, [1 0 0]);
        axis square;
        plot(xlim, [0 0], '--', 'color', [.5 .5 .5]);
        xlabel('From response (s)');
        set(gca, 'YAxisLocation', 'right');
        ylim2 = get(gca,'ylim');
    end
    ylim_match = [min([ylim1(1),ylim2(1)]), max([ylim1(2), ylim2(2)])];
    ax = get(gcf,'Children');
    for n = 1:length(ax)
        set(ax(n), 'ylim', ylim_match);
    end
else
    fh = [];
end

stat.t = t;
stat.stkernel1 = stkernel1;
stat.stkernel2 = stkernel2;
stat.rekernel1 = rekernel1;
stat.rekernel2 = rekernel2;

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
lh = plot(x,y,'Color', color, 'LineWidth', .5);

end


function [coh_task1, coh_task2, coh_main, coh_orth] = get_coh(morph_level, cond)

% coh for task 1/2
morph_level = cell2mat(morph_level);
coh_task1 = morph_level(:,1);
coh_task2 = morph_level(:,2);

% coh main/orth
coh_main = nan(length(coh_task1),1);
coh_orth = nan(length(coh_task1),1);
coh_main(cond==1) = coh_task1(cond==1); coh_orth(cond==1) = coh_task2(cond==1); 
coh_main(cond==2) = coh_task2(cond==2); coh_orth(cond==2) = coh_task1(cond==2);

end


function [kernel_mean, kernel_sem] = calKernelMean(kernel)

nsubj = length(kernel);
nfeature = length(kernel{1});
medRT = length(kernel{1}{1});

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


function C = nanconv(A,B,shape)

if (nargin < 3)
  shape = 'full';
end

% N = convn(A,B,shape);
% D = convn(ones(size(A)),B,shape);
% C = N./D;

AO = ones(size(A));
I = isnan(A);
A(I) = 0;
AO(I) = 0;
C = convn(A,B,shape) ./ convn(AO,B,shape);
C(I) = nan;

end


function se = calc_diff_SE(data1, data2)

se1 = nanstd(data1,0,1) ./ sqrt(sum(~isnan(data1),1));
se2 = nanstd(data2,0,1) ./ sqrt(sum(~isnan(data2),1));
se = sqrt(se1.^2 + se2.^2);

end