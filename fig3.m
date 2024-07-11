clear;
% opengl software
addpath('./functions/');

D = load('./data/preproc/context_dependent_face_categorization_task.mat').trial_data;



%% Fig 3A
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
morph_level = cellfun(@(x) x.morph_level, D, 'uni', 0);
cond = cellfun(@(x) x.cond, D);
fluc = cellfun(@(x) x.fluc, D, 'uni', 0);
resp = cellfun(@(x) x.resp, D);

clear opt;
opt.smoothing = 3;
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.fit_file = cellfun(@(x) ['./data/fit/nonswitch/', x, '_DDMfit.mat'], opt.subj_list, 'uni', 0);
opt.ylim = [-1 5.5];
run_show_kernel(cond_switch, morph_level, cond, fluc, resp, subj, opt);

%% Fig 3B
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
morph_level = cellfun(@(x) x.morph_level, D, 'uni', 0);
cond = cellfun(@(x) x.cond, D);
fluc = cellfun(@(x) x.fluc, D, 'uni', 0);
resp = cellfun(@(x) x.resp, D);

clear opt;
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.fit_file = cellfun(@(x) ['./data/fit/nonswitch/', x, '_DDMfit.mat'], opt.subj_list, 'uni', 0);
opt.ylim = [-1 5.5];
run_show_kernel_orthCoh(cond_switch, morph_level, cond, fluc, resp, subj, opt);

%% Fig 3D
clear opt;
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.pos_cho = 2;
opt.ylim = [.6 1.6];

D_ns = preformat_fitresult_for_plotting('./data/fit', 'nonswitch', '_DDMfit.mat', opt);
coh_ns = D_ns.coh;
resp_ns = D_ns.resp;
fit_pchoice1_ns = D_ns.fit_pchoice1;
rt_ns = D_ns.rt;
fit_rt_ns = D_ns.fit_rt;
subj_ns = D_ns.subj;

D_s = preformat_fitresult_for_plotting('./data/fit', 'switch', '_DDMfit.mat', opt);
coh_s = D_s.coh;
resp_s = D_s.resp;
fit_pchoice1_s = D_s.fit_pchoice1;
rt_s = D_s.rt;
fit_rt_s = D_s.fit_rt;
subj_s = D_s.subj;
opt.param_final = D_s.param_final;

run_show_fit_choiceRT_comparison(coh_ns, resp_ns, fit_pchoice1_ns, rt_ns, fit_rt_ns, subj_ns, ...
    coh_s, resp_s, fit_pchoice1_s, rt_s, fit_rt_s, subj_s, ...
    opt);

%% Fig 3E
sim_file = './data/fit/switch/sim/rt_CDF.mat';
run_show_rtCDF_allCoh_comparison(sim_file);

%% Fig 3F
sim_file = './data/fit/nonswitch/sim/kernel.mat';

clear opt;
opt.filter_cond_switch = 1;
opt.ylim = [-1 5.5];
run_show_fit_kernel(sim_file, opt);

%% Fig 3G
sim_file = './data/fit/switch/sim/kernel.mat';

clear opt;
opt.filter_cond_switch = 2;
opt.ylim = [-1 5.5];
run_show_fit_kernel(sim_file, opt);

