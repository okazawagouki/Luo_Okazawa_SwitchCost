function [fh, stat] = run_show_rtCDF_allCoh_comparison(sim_file)

stat = load(sim_file).stat;
fh = show_rtCDF(stat);

end



function [fh, stat] = show_rtCDF(stat)

opt.prctile = [3 10 20 30 40 50 60 70 80 90];
opt.coh_range = [0    0.12;
                 0.24 0.96];

%% data
data_rt_cdf_ns = stat.data_rt_cdf_ns;
data_rt_cdf_s = stat.data_rt_cdf_s;

%% model
model_rt_ns = stat.model_rt_ns;
model_rt_s = stat.model_rt_s;

%% plot
fh = figure('color','w','Position',[100 100 350 150]);
for n = 1:size(opt.coh_range,1)
    subplot(1,2,n);
    hold on
    m1 = cdfplot(model_rt_ns{n});
    set(m1,'color', .5+.5*[0 0 0], 'linew', 1.5);
    d1 = plot(data_rt_cdf_ns{n}, opt.prctile/100, '.k', 'MarkerSize', 7);
    m2 = cdfplot(model_rt_s{n});
    set(m2,'color', .5+.5*[1 0 0], 'linew', 1.5);
    plot(data_rt_cdf_s{n}, opt.prctile/100, '.r', 'MarkerSize', 7);
    grid off
    title([num2str(opt.coh_range(n,1)*100,'%.0f'), '-', num2str(opt.coh_range(n,2)*100,'%.0f'), '% Morph'])
    xlabel('Reaction time (s)')
    ylabel({'Cumulative','distribution'})
    xlim([0 3])
    if n==1
        legend([d1 m1], {'Data','Model'}, 'Location', 'northwest'); legend boxoff;
    end
end

%% save
stat.data_rt_cdf_ns = data_rt_cdf_ns;
stat.model_rt_ns = model_rt_ns;
stat.data_rt_cdf_s = data_rt_cdf_s;
stat.model_rt_s = model_rt_s;

end