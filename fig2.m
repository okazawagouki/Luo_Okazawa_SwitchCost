clear;
% opengl software
addpath('./functions/');

D = load('./data/preproc/context_dependent_face_categorization_task.mat').trial_data;



%% Fig 2A: Psychometric and chronometric functions along task axis
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);
rt = cellfun(@(x) x.rt, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch', 'Switch'};
opt.xlabel = {'Task axis', '(% Morph)'};
opt.ylim = [0.6, 1.6];
run_show_choiceRT(cond_switch, coh, resp, rt, subj, opt);

%% Fig 2B: Psychometric and chronometric functions along orthogonal axis
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh_orth = cellfun(@(x) x.coh_orth, D);
resp = cellfun(@(x) x.resp, D);
rt = cellfun(@(x) x.rt, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch', 'Switch'};
opt.xlabel = {'Orthogonal axis', '(% Morph)'};
opt.ylim = [0.6, 1.6];
run_show_choiceRT(cond_switch, coh_orth, resp, rt, subj, opt);

%% Fig 2C
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
rt = cellfun(@(x) x.rt, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.xlabel = {'Task axis', '(Absolute % Morph)'};
opt.ylim = [0 0.4];
run_show_drtCoh_unsigned(cond_switch, coh, rt, subj, opt);

%% Fig 2D
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh_orth = cellfun(@(x) x.coh_orth, D);
rt = cellfun(@(x) x.rt, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.xlabel = {'Orthogonal axis', '(Absolute % Morph)'};
opt.ylim = [0 0.4];
run_show_drtCoh_unsigned(cond_switch, coh_orth, rt, subj, opt);

%% Fig 2E
subj = cellfun(@(x) x.subj, D, 'uni', 0);
n_after_switch = cellfun(@(x)x.n_after_switch, D);
cond_switch = cellfun(@(x) x.cond_switch, D);
rt = cellfun(@(x) x.rt, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.ylim = [-0.08 0.25];
opt.xlim = [-.5 5.5]+1;
run_show_drtAfterNSwitch(cond_switch, n_after_switch, rt, subj, opt);

%% Fig 2F
subj = cellfun(@(x) x.subj, D, 'uni', 0);
n_after_nonswitch = cellfun(@(x)x.n_after_nonswitch, D);
cond_switch = cellfun(@(x) x.cond_switch, D);
rt = cellfun(@(x) x.rt, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.ylim = [-0.08 0.4];
opt.xlim = [.5 4.5];
run_show_drtAfterNNonswitch(cond_switch, n_after_nonswitch, rt, subj, opt);

%% Fig 2G
csi = cellfun(@(x) x.csi, D);

fh = figure('color','w','Position',[100 100 350 150]);
subplot(1,2,1);
hold on;
histogram(csi, 15, 'Normalization', 'probability',...
    'EdgeColor', 'k', 'EdgeAlpha', 1, 'FaceColor', 'none', 'FaceAlpha', 0);
xlabel('CSI (s)')
ylabel('Probability')

%% Fig 2H (left panel)
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);
rt = cellfun(@(x) x.rt, D);
csi = cellfun(@(x) x.csi, D);
I = csi<=0.72;

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch', 'Switch'};
opt.xlabel = {'Task axis', '(% Morph)'};
opt.ylim = [0.6, 1.6];
run_show_choiceRT(cond_switch(I), coh(I), resp(I), rt(I), subj(I), opt);

%% Fig 2H (right panel)
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);
rt = cellfun(@(x) x.rt, D);
csi = cellfun(@(x) x.csi, D);
I = csi>0.72;

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch', 'Switch'};
opt.xlabel = {'Task axis', '(% Morph)'};
opt.ylim = [0.6, 1.6];
run_show_choiceRT(cond_switch(I), coh(I), resp(I), rt(I), subj(I), opt);

%% Fig 2I
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
rt = cellfun(@(x) x.rt, D);
csi = cellfun(@(x) x.csi, D);
coh = cellfun(@(x) x.coh, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.ylim = [0 0.4];
opt.nbin = 7;
opt.plot_idv = 0;
run_show_drt_csi(cond_switch, coh, csi, rt, subj, opt);

%% Fig 2J (upper panels)
subj = cellfun(@(x) x.subj, D, 'uni', 0);
morph_flip = cellfun(@(x) x.morph_flip, D, 'uni', 0);
cond = cellfun(@(x) x.cond, D);
resp2_nonflip = cellfun(@(x) x.resp_nonflip==2, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
run_show_heatmap(resp2_nonflip, morph_flip, cond, subj, opt);

%% Fig 2J (lower panels)
subj = cellfun(@(x) x.subj, D, 'uni', 0);
morph_flip = cellfun(@(x) x.morph_flip, D, 'uni', 0);
cond = cellfun(@(x) x.cond, D);
rt = cellfun(@(x) x.rt, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
run_show_heatmap(rt, morph_flip, cond, subj, opt);




