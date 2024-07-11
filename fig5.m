clear;
% opengl software
addpath('./functions/');

D = load('./data/preproc/fixed_stimulus_duration_task.mat').trial_data;



%% Fig 5B, D
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);
targ_cor = cellfun(@(x)x.targ_cor, D);
rt = cellfun(@(x) x.rt, D);

clear opt
opt.subj_list = {'005','012','013','038','020','057','054'};
run_show_choiceRT_unsigned(cond_switch, coh, resp, targ_cor, rt, subj, opt);

%% Fig 5C
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x)x.cond_switch, D);
coh = cellfun(@(x)x.coh, D);
resp = cellfun(@(x)x.resp, D);
targ_cor = cellfun(@(x)x.targ_cor, D);
csi = cellfun(@(x)x.csi, D) * 1e3; % s to ms

clear opt
opt.subj_list = {'005','012','013','038','020','057','054'};
opt.stimdur_or_csi = 'csi';
opt.nbin = 4;
run_show_dthreshold_logistic(cond_switch, coh, resp, targ_cor, csi, subj, opt);

