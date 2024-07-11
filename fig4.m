clear;
% opengl software
addpath('./functions/');

D = load('./data/preproc/context_dependent_face_categorization_task.mat').trial_data;



%% Fig 4B
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

D_s = preformat_fitresult_for_plotting('./data/fit', 'T0', '_DDMfit.mat', opt);
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

%% Fig 4C
sim_file = './data/fit/T0/sim/kernel.mat';

clear opt;
opt.filter_cond_switch = 2;
opt.ylim = [-1 5.5];
run_show_fit_kernel(sim_file, opt);

%% Fig 4E (left panel)
sim_file = './data/fit/prep/sim/kernel_shortCSI.mat';

clear opt;
opt.filter_cond_switch = 2;
opt.ylim = [-1 5.5];
run_show_fit_kernel(sim_file, opt);

%% Fig 4E (right panel)
sim_file = './data/fit/prep/sim/kernel_longCSI.mat';

clear opt;
opt.filter_cond_switch = 2;
opt.ylim = [-1 5.5];
run_show_fit_kernel(sim_file, opt);

%% Fig 4F
clear opt
opt.ylim = [-0.05, 0.3];

sim_file = './data/fit/prep/sim/stat_drt_csi.mat';
run_show_fit_drt_csi(sim_file, opt);

%% Fig 4G
ProjectDirCompare = {'switch','T0', 'prep', 'leak','B','k','Bk'};
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.FormalProjectName = {'Non-decision time\uparrow', 'Task preparation', 'Leak', 'Bound\uparrow', 'Sensitivity\downarrow', 'Bound\uparrow & sensitivity\downarrow'};
get_summary_LL(ProjectDirCompare, './data/fit', opt);


