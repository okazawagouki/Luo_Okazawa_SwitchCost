function D = preformat_fitresult_for_plotting(MainInterimDir, ProjectDir, postfix, opt)

% example:
% D = preformat_fitresult_for_plotting(MainInterimDir, ProjectDir, '_DDMfluc_fit_largeDt.mat', opt)

for s = 1:length(opt.subj_list)
    bestfitfile{s} = fullfile(MainInterimDir, ProjectDir, [opt.subj_list{s}, postfix]);
    fitresult{s} = load(bestfitfile{s}).fitresult;
    param_final{s} = fitresult{s}.model_param.final;
end

coh = cellfun(@(x) x.trial_data.coh, fitresult, 'uni', 0);                   coh = coh(:);                   coh = cell2mat(coh);
resp = cellfun(@(x) x.trial_data.response, fitresult, 'uni', 0);             resp = resp(:);                 resp = cell2mat(resp);
targ_cor = cellfun(@(x) x.trial_data.targ_cor, fitresult, 'uni', 0);         targ_cor = targ_cor(:);         targ_cor = cell2mat(targ_cor);
rt = cellfun(@(x) x.trial_data.rt, fitresult, 'uni', 0);                     rt = rt(:);                     rt = cell2mat(rt);
fit_pchoice1 = cellfun(@(x) x.trial_data.fit_pchoice1, fitresult, 'uni', 0); fit_pchoice1 = fit_pchoice1(:); fit_pchoice1 = cell2mat(fit_pchoice1);
fit_rt = cellfun(@(x) x.trial_data.fit_rt, fitresult, 'uni', 0);             fit_rt = fit_rt(:);             fit_rt = cell2mat(fit_rt);
subj = [];
for s = 1:length(opt.subj_list)
    subj = [subj; repmat(opt.subj_list(s), size(fitresult{s}.trial_data.coh))];
end

D = struct('coh', coh, 'resp', resp, 'targ_cor', targ_cor, 'rt', rt, ...
    'fit_pchoice1', fit_pchoice1, 'fit_rt', fit_rt);
D.subj = subj;
D.param_final = param_final;

end