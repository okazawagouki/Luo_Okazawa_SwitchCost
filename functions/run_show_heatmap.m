function [fh, fh_idv] = run_show_heatmap(respRT, morph_flip, cond, subj, opt)

opt.plot_idv_fh = false;

nsubj = length(opt.subj_list);
for s = 1:nsubj
    I = strcmp(subj, opt.subj_list{s});
    fh_idv{s} = show_heatmap(cond(I), morph_flip(I), respRT(I), opt);
end

opt.plot_idv_fh = true;
fh = show_heatmap(cond, morph_flip, respRT, opt);

end



function fh = show_heatmap(cond, morph_flip, respRT, opt)

% averaging
morph_task1 = cellfun(@(x) x(1), morph_flip);
morph_task2 = cellfun(@(x) x(2), morph_flip);

mean_respRT{1} = calMean2Dim(respRT(cond==1), morph_task1(cond==1), morph_task2(cond==1));
mean_respRT{2} = calMean2Dim(respRT(cond==2), morph_task1(cond==2), morph_task2(cond==2));

% plot
if opt.plot_idv_fh
fh = figure('color','w','Position',[100 100 350 150]);
    for n = 1:2
        ax(n) = subplot(1,2,n);
        imagesc(mean_respRT{n});
        if isequal(unique(respRT), [0;1])
            colormap(winter)
        else
            colormap(spring)
        end
        
        axis equal; axis off; axis on
        xlabel({'Context 1 axis', '(% Morph)'}); ylabel({'Context 2 axis', '(% Morph)'})
        title(['Context ', num2str(n)])
        xlim([0.5 11.5])
        ylim([0.5 11.5])
        XAxis.Visible = 'off'; YAxis.Visible = 'off'; box off
        xticks([]); yticks([])
    end
    colorbar;
    set(ax(1),'Position',[0.1 0.3 0.3 0.6],'XColor','k','YColor','k')
    set(ax(2),'Position',[0.6 0.3 0.3 0.6],'XColor','k','YColor','k')
else
    fh = [];
end

end



function mean_data = calMean2Dim(data, cond1, cond2)

%MEAN_RESP2 group the mean probability of choosing target 2, based on morph
%level of two dimensions

% input ==>
% resp2_nonflip ..0: does not choose target 2; 1: choose target 2
% morph_task1 ..morph level of task 1
% morph_task2 ..morph level of task 2

% output ==>
% mean_resp2 ..size=[ncoh2, ncoh1], x-axis: morph_task1 from -100 to 100,
% y-axis: morph_task2 from -100 to 100, origin: lower left

ncond1 = length(unique(cond1));
ncond2 = length(unique(cond2));

tb = table(data, cond1, cond2);
tbstat = grpstats(tb, ["cond1", "cond2"]);
mean_data = tbstat.mean_data;
mean_data = reshape(mean_data, [ncond2, ncond1]); % reshape(mean_resp2, [row column]);
mean_data = flip(mean_data); % flip according to row

end