clear;
addpath('./functions/');



%% Psychometric and chronometric function
clear;

S = load('./data/preproc/data.mat');
subj = S.subj;
cond_switch = S.cond_switch;
cond = S.cond;
morph_level = S.morph_level;
coh = nan(size(subj));
for n = 1:length(subj)
    coh(n) = morph_level{n}(cond(n));
end
resp = S.resp;
rt = S.rt;
I = ~isnan(rt);

clear opt
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.color = {'k', 'r'};
opt.legend = {'Non-switch trial', 'Switch trial'};
opt.xlabel = {'% Morph along', 'task axis'};
opt.ylim = [0.6, 1.6];
run_show_choiceRT(cond_switch(I), coh(I), resp(I), rt(I), subj(I), opt);



%% Psychophysical kernel
clear;

S = load('./data/preproc/data.mat');
subj = S.subj;
cond_switch = S.cond_switch;
morph_level = S.morph_level;
cond = S.cond;
fluc = S.fluc;
resp = S.resp;
rt = S.rt;
I = ~isnan(rt);

clear opt;
opt.subj_list = {'002','005','007','008','012','013','010','014'};
opt.fit_file = cellfun(@(x) ['./data/fit/', x, '_DDMfluc_fit_smallDt.mat'], opt.subj_list, 'uni', 0);
opt.ylim = [-1 5.5];
run_show_kernel(cond_switch(I), morph_level(I), cond(I), fluc(I), resp(I), subj(I), opt);



