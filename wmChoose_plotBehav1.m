% wmChoose_plotBehav1.m
%
% plot general results from wmChoose behavioral experiment


function wmChoose_plotBehav1(subj)

root = '/Volumes/data/wmChoose/';
if nargin < 1 || isempty(subj)
    %subj = {'KD','CC','EK','MR','AB'};
    subj = {'aa','ab1','ac1','ac2','ae','af','ag'}; %aa1

end


f_err_thresh = 5;

% load everything
% [[ for now, just err for i_sacc, f_sacc; we'll augment this later on ]]

all_conds = [];
all_i_err = [];
all_f_err = [];
all_ct = []; % chosen target

all_subj = [];

for ss = 1:length(subj)
    
    fn = sprintf('%s/data/%s_wmChoose_behav.mat',root,subj{ss});
    thisdata = load(fn);
    
    all_conds = [all_conds;thisdata.c_all(:,1)];
    all_ct = [all_ct;thisdata.s_all.ti_all(:,6)]; % chosen target
    
    all_i_err = [all_i_err;thisdata.s_all.i_sacc_err];
    all_f_err = [all_f_err;thisdata.s_all.f_sacc_err];
    
    all_subj = [all_subj;ss*ones(size(thisdata.c_all,1),1)];

    
    clear thisdata;
    
end


% first, dumbest possible figure: for each subj, conds 1 vs 2 vs 3
figure;


% save mean across trials within each subj
this_m_i = nan(length(subj),3);
this_m_f = nan(length(subj),3);
for ss = 1:length(subj)
    
    for cc = 1:3
        
        % TODO: and error is < 5 or something...
        thisidx = all_subj==ss & all_conds==cc & all_f_err<f_err_thresh;
        this_m_i(ss,cc) = nanmean(all_i_err(thisidx));
        this_m_f(ss,cc) = nanmean(all_f_err(thisidx));
        
        clear thisidx;
        
    end
end

subplot(1,2,1); hold on;
for ss = 1:length(subj)
    plot(1:3,this_m_i(ss,:),'-','Color',[0.5 0.5 0.5]);
end
plot(1:3,mean(this_m_i,1),'k-','LineWidth',2);
title('Initial saccade');
set(gca,'XTick',[1:3],'XTickLabel',{'R1','R2','choose'});

subplot(1,2,2); hold on;
for ss = 1:length(subj)
    plot(1:3,this_m_f(ss,:),'-','Color',[0.5 0.5 0.5]);
end
plot(1:3,mean(this_m_f,1),'k-','LineWidth',2);
title('Final saccade');
set(gca,'XTick',[1:3],'XTickLabel',{'R1','R2','choose'});



return