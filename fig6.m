clear;
% opengl software
addpath('./functions/');

D_main = load('./data/preproc/context_dependent_face_categorization_task.mat').trial_data;
D_idt_col = load('./data/preproc/identity_versus_color.mat').trial_data;
D_mot_col = load('./data/preproc/motion_versus_color.mat').trial_data;



%% Fig 6D (left panels, same as Fig 2A)
D = D_main;
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);
rt = cellfun(@(x) x.rt, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch', 'Switch'};
opt.xlabel = 'Stimulus strength (%)';
opt.ylim = [0.5, 1.7];
run_show_choiceRT(cond_switch, coh, resp, rt, subj, opt);


%% Fig 6D (middle panels)
D = D_idt_col;
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);
rt = cellfun(@(x) x.rt, D);

clear opt
opt.subj_list = {'020','048','049','037','038','057','012'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch', 'Switch'};
opt.xlabel = 'Stimulus strength (%)';
opt.ylim = [0.5, 1.7];
run_show_choiceRT(cond_switch, coh, resp, rt, subj, opt);

%% Fig 6D (right panels)
D = D_mot_col;
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);
rt = cellfun(@(x) x.rt, D);
cond = cellfun(@(x) x.cond, D);
coh = coh_normalization(coh, cond); % coherence normalization

clear opt
opt.subj_list = {'025','027','017','053','050','014','054'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch', 'Switch'};
opt.xlabel = 'Relative stimulus strength';
opt.xlim = [-1.1 1.1]; % relative stimulus strength
opt.ylim = [0.5, 1.7];
run_show_choiceRT_mot_col(cond_switch, coh, resp, rt, subj, opt);

%% Fig 6E (left panels)
D = D_main;
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);
rt = cellfun(@(x) x.rt, D);
cond = cellfun(@(x) x.cond, D);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch', 'Switch'};
opt.xlabel = 'Stimulus strength (%)';
opt.ylim = [0.5, 1.7];

I = cond==1;
run_show_choiceRT(cond_switch(I), coh(I), resp(I), rt(I), subj(I), opt); % Identity trials

I = cond==2;
run_show_choiceRT(cond_switch(I), coh(I), resp(I), rt(I), subj(I), opt); % Expression trials

%% Fig 6E (middle panels)
D = D_idt_col;
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);
rt = cellfun(@(x) x.rt, D);
cond = cellfun(@(x) x.cond, D);

clear opt
opt.subj_list = {'020','048','049','037','038','057','012'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch', 'Switch'};
opt.xlabel = 'Stimulus strength (%)';
opt.ylim = [0.5, 1.7];

I = cond==1;
run_show_choiceRT(cond_switch(I), coh(I), resp(I), rt(I), subj(I), opt); % Identity trials

I = cond==2;
run_show_choiceRT(cond_switch(I), coh(I), resp(I), rt(I), subj(I), opt); % Color trials

%% Fig 6E (right panels)
D = D_mot_col;
subj = cellfun(@(x) x.subj, D, 'uni', 0);
cond_switch = cellfun(@(x) x.cond_switch, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);
rt = cellfun(@(x) x.rt, D);
cond = cellfun(@(x) x.cond, D);
coh = coh_normalization(coh, cond); % coherence normalization

clear opt
opt.subj_list = {'025','027','017','053','050','014','054'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch', 'Switch'};
opt.xlabel = 'Relative stimulus strength';
opt.xlim = [-1.1 1.1]; % relative stimulus strength
opt.ylim = [0.5, 1.7];

I = cond==1;
run_show_choiceRT_mot_col(cond_switch(I), coh(I), resp(I), rt(I), subj(I), opt); % Motion trials

I = cond==2;
run_show_choiceRT_mot_col(cond_switch(I), coh(I), resp(I), rt(I), subj(I), opt); % Color trials








