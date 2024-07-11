function fh = get_summary_LL(ProjectDirCompare, MainInterimDir, opt)

% load fitting result
nmodel = length(ProjectDirCompare);
nsubj = length(opt.subj_list);
LL = nan(nmodel, nsubj);
nparam = nan(nmodel, nsubj);
ndata = nan(nmodel, nsubj);
for m = 1:length(ProjectDirCompare)
    for s = 1:length(opt.subj_list)
        file_path = fullfile(MainInterimDir, ProjectDirCompare{m}, [opt.subj_list{s}, '_DDMfit.mat']);
        S = load(file_path).fitresult;
        LL(m,s) = S.modelLL;
        nparam(m,s) = sum(S.model_param.fixed==0);
        ndata(m,s) = length(S.trial_data.LL);
    end
end

% summary
BIC = log(ndata).*nparam - 2*LL; BIC_mean = mean(BIC,2); % BIC = ln(ndata)*nparam - 2*ln(L)
dBIC = BIC_mean(2:end) - BIC_mean(1);

% plot dBIC
fh = figure('Color', 'w', 'Position', [200 200 220 200]);
subplot(1,2,1)
bar(1:nmodel-1, dBIC, .8, 'FaceColor', '#A9A9A9')
axis normal
box off
ylim([0 120])
yticks(linspace(0,120,5))
ylabel('\DeltaBIC')
xticks(1:nmodel-1)
xticklabels(opt.FormalProjectName)
xtickangle(70)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',7)

end
