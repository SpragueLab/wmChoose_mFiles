% wmChoose_plotResults.m
%
% general plotting script for wmChoose - plots sacc metrics as a function
% of condition (R1, R2-cued, R2-choose)
%
% TCS 4/13/2018


root = '/Volumes/data/wmChoose/';

subj = {'aa1','aa2','ab1','ab2','ac1','ac2','ae','af','ag'};

WHICH_EXCL = [11 13 20 21 22]; % don't exclude trials w/ calibration failures for now...

% for now, let's use cat_struct to load/concatenate all data...
all_subj = nan(1000*length(subj),1);
u_subj = unique(cellfun(@(s) s(1:2),subj,'uniformoutput',0));

TARG_ECC = 12;

niter = 1000;

all_data = [];

startidx = 1;

for ss = 1:length(subj)
    fn = sprintf('%s/data/%s_wmChoose_behav.mat',root,subj{ss});
    fprintf('Loading %s\n',fn);
    this_data = load(fn);
    
    this_subj = find(strcmpi(u_subj,subj{ss}(1:2)));
    
    all_data = cat_struct(all_data,this_data);
    all_subj(startidx:(startidx-1+size(this_data.c_all,1))) = this_subj;
    
    startidx = startidx+size(this_data.c_all,1);
    
    clear this_subj this_data;
end

% let's try this pattern for now
all_subj = all_subj(1:(startidx-1));
all_data.subj_all = all_subj;

% determine which trials to include
% first, narrow based on saccade preprocessing/scoring exclusions
% (wmChoose_extractSaccadeData1.m)
all_data.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data.s_all.excl_trial, 'UniformOutput',false));

% drop trials with very short (< 100 ms) or very long RT (> 1 s)
all_data.use_trial(all_data.s_all.i_sacc_rt<0.1 | all_data.s_all.i_sacc_rt>1.5) = 0;

%% first, plot mean i_sacc, f_sacc error as a function of condition

mean_fig = figure;
scatter_fig = figure;

to_plot = {'i_sacc_err','f_sacc_err','i_sacc_rt'};
cu = unique(all_data.c_all(:,1));

cond_str = {'R1','R2-cue','R2-choose'};

cond_colors = lines(length(cu));

cond_pairs = [1 2; 2 3; 1 3]; % x, y axes of scatterplot


for pp = 1:length(to_plot)
    figure(mean_fig);
    subplot(1,length(to_plot),pp); hold on;
    
    thisd = nan(length(u_subj),length(cu));
    for cc = 1:length(cu)
        for ss = 1:length(u_subj)
            thisidx = all_data.subj_all==ss & all_data.c_all(:,1)==cu(cc) & all_data.use_trial==1;
            thisd(ss,cc) = mean(all_data.s_all.(to_plot{pp})(thisidx));
        end
        plot(cc+[-0.35 0.35],[1 1]*mean(thisd(:,cc)),'-','LineWidth',2.5,'Color',cond_colors(cc,:))
    end
    plot(1:length(cu),thisd.','-','Color',[0.5 0.5 0.5]);
    
    set(gca,'XTick',1:length(cu),'TickDir','out','LineWidth',1.5,'XTickLabel',cond_str,'FontSize',14,'XTickLabelRotation',-45);
    xlim([0.5 0.5+length(cu)]);
    title(to_plot{pp},'Interpreter','none');
    
    figure(scatter_fig);
    for cp = 1:size(cond_pairs,1)
        subplot(size(cond_pairs,1),length(to_plot),pp+(cp-1)*length(to_plot)); hold on;
        plot(thisd(:,cond_pairs(cp,1)),thisd(:,cond_pairs(cp,2)),'o','LineWidth',1.5,'MarkerSize',5,'Color',[0.5 0.5 0.5]);
        plot([0 3],[0 3],'k--','LineWidth',1.5);
        if pp==3 % if RT use different xlim,ylim
            xlim([0 1]); ylim([0 1]);
        else
            xlim([0 3]); ylim([0 3]);
        end
        axis square;
        xlabel(cond_str{cond_pairs(cp,1)},'Interpreter','none');
        ylabel(cond_str{cond_pairs(cp,2)},'Interpreter','none');
        if cp == 1
            title(to_plot{pp},'Interpreter','none');
        end
        set(gca,'TickDir','out','FontSize',14,'LineWidth',1.5);
    end
    
end

%% 2d distribution of all trials for i_sacc, f_sacc for each subj, condition
to_plot_2d = {'i_sacc','f_sacc'};

dist_2d_figs = nan(length(to_plot_2d),1);

for pp = 1:length(to_plot_2d)
    dist_2d_figs(pp) = figure;
    for cc = 1:length(cu)
        for ss = 1:length(u_subj)
            
            subplot(length(cu),length(u_subj),ss+(cc-1)*length(u_subj)); hold on;
            
            thisidx = all_data.subj_all==ss & all_data.c_all(:,1)==cu(cc) & all_data.use_trial==1;
            %plot(all_data.s_all.(to_plot_2d{pp})(thisidx,1),all_data.s_all.(to_plot_2d{pp})(thisidx,2),'.','Color',cond_colors(cc,:));
            plot(0,0,'ko','MarkerFaceColor','k','MarkerSize',6);   % TODO: replace w/ scatter, 'MarkerFaceAlpha' = 0.5
            scatter(all_data.s_all.(to_plot_2d{pp})(thisidx,1),all_data.s_all.(to_plot_2d{pp})(thisidx,2),20,cond_colors(cc,:),'filled','MarkerFaceAlpha',0.2);
            plot(TARG_ECC,0,'k+','LineWidth',1.5,'MarkerSize',5);
            plot([0 3],[-5 -5],'k-','LineWidth',1.5);
            axis equal off;
            
            % if first subj, draw 'YLabel'
            if ss == 1
                text(-15,0,cond_str{cc},'FontSize',12,'HorizontalAlignment','center','Rotation',90);
            end
            
            % if first condition, draw title
            if cc == 1
                %title(u_subj{ss});
                text(5,12,u_subj{ss},'HorizontalAlignment','center','FontSize',12);
            end
        end
    end
    
    axes('Position',[0.45 0.95 0.1 0.05]);
    text(0.5,0.5,to_plot_2d{pp},'HorizontalAlignment','center','FontSize',14,'FontWeight','bold','Interpreter','none');
    axis off;
    %match_xlim(get(gcf,'Children'));
    %match_ylim(get(gcf,'Children'));
    set(get(gcf,'Children'),'XLim',[-5 15],'YLim',[-15 15]);axis equal;
end



%% std dev (radial, tangential) and distributions for each saccade (i_sacc, f_sacc), condition

% use to_plot_2d again
% 4 histogram figures:
% - response x rad/tang, each:
% - each row is a condition, each col a subj
% - each cell a histogram of radial/tangential i_sacc/f_sacc error
%   distribution

% and also a summary plot (initial, final, radial, tangential std dev by
% condition)

dim_to_plot = [1 2]; % first or second dimension of s_all.i_sacc, f_sacc
dim_targ = [TARG_ECC 0]; % what to subtract from i_sacc, f_sacc for each dim
dim_str = {'Radial','Tangential'}; % x, y after rotating

this_nrows = length(dim_to_plot)*length(to_plot_2d);

all_std = nan(length(dim_to_plot),length(to_plot_2d),length(cu),length(u_subj)); % store std dev for each condition; subj (dumb way for now...)

fig_dist = figure;
for pp = 1:length(to_plot_2d)
    for dd = 1:length(dim_to_plot)
        
        thisrow = (pp-1)*length(dim_to_plot)+dd;

        %figs_dist(pp,dd) = figure;
        
        for ss = 1:length(u_subj)
            subplot(this_nrows,length(u_subj),ss+(thisrow-1)*length(u_subj)); hold on;

            for cc = 1:length(cu)
                
                
                %subplot(length(cu),length(u_subj),ss+(cc-1)*length(u_subj)); hold on;
                
                
                thisidx = all_data.subj_all==ss & all_data.c_all(:,1)==cu(cc) & all_data.use_trial==1;
                
                [thish,thise] = histcounts(all_data.s_all.(to_plot_2d{pp})(thisidx,dim_to_plot(dd)) - dim_targ(dd),'BinWidth',0.5,'Normalization','pdf');
                plot(mean([thise(1:end-1);thise(2:end)],1),thish,'-','LineWidth',1,'Color',cond_colors(cc,:));
                
                all_std(dd,pp,cc,ss) = std(all_data.s_all.(to_plot_2d{pp})(thisidx,dim_to_plot(dd)) - dim_targ(dd));
                
                clear thish thise;
                
                if ss == 1
                    ylabel(sprintf('%s - %s',dim_str{dd},to_plot_2d{pp}),'Interpreter','none');
                end
                
                if thisrow==1
                    title(u_subj{ss});
                end
                
                clear thisidx;
                xlim([-7 7]);
                
            end
        end
        
        
    end
end

match_xlim(get(gcf,'Children'));
match_ylim(get(gcf,'Children'));

set(get(gcf,'Children')','YTickLabel',[]);


% now plot all the std devs
fig_std_sum = figure;
for pp = 1:length(to_plot_2d)
    for dd = 1:length(dim_to_plot)
        
        subplot(length(to_plot_2d),length(dim_to_plot),dd+(pp-1)*length(to_plot_2d)); hold on;
        
        plot(1:length(cu),squeeze(all_std(dd,pp,:,:)),'o-','Color',[0.5 0.5 0.5]);
        
        for cc = 1:length(cu)
            
            thismu = mean(squeeze(all_std(dd,pp,cc,:)));
            
            plot(cc+[-0.35 0.35],thismu*[1 1],'-','LineWidth',2,'Color',cond_colors(cc,:));
            
            clear thismu;
        end
        
        set(gca,'XTick',1:length(cu),'XTickLabel',cond_str,'XTickLabelRotation',-45);
        xlim([0.35 3.65]);
        
        if dd == 1
            ylabel(to_plot_2d{pp},'Interpreter','none');
        end
        
        if pp == 1
            title(dim_str{dd});
        end
        
    end
end
match_ylim(get(fig_std_sum,'Children'));

% TODO: 
% - trial exclusion % by condition (maybe broken down by exclusion type?)
% - trace for each condition; distribution for each condition

%% stats - shuffle condition labels within each subj before computing distributions, F-scores
% - use only included trials? yes

% set the random number generator before we do stats
rng(wmChoose_randSeed());

allF = nan(length(to_plot_2d),length(dim_to_plot),niter+1);   % params x dimensions x iterations?
allT = cell(size(cond_pairs,1),1);  % each of these same dims as above
for cp_idx = 1:size(cond_pairs,1)
    allT{cp_idx} = nan(length(to_plot_2d),length(dim_to_plot),niter+1);
end

all_labels = [all_data.c_all(all_data.use_trial==1,1) all_data.subj_all(all_data.use_trial==1)]; % 1 column for condition, 1 for subj

% first iteration is 'real' (only shuffle on 2nd-nth iteration

for ii = 1:(niter+1)
    
    % if ii ~= 1, shuffle labels within each subj
    if ii~=1
        this_labels = nan(size(all_labels));
        this_labels(:,2) = all_labels(:,2);
        for ss = 1:length(u_subj)
            tmplabel = all_labels(all_labels(:,2)==ss,1);
            this_labels(all_labels(:,2)==ss,1) = tmplabel(randperm(length(tmplabel)));
        end
    else
        this_labels = all_labels;
    end
    
    
    % for saccade parameters
    for pp = 1:length(to_plot_2d)
        
        % radial, tangential error
        for dd = 1:length(dim_to_plot)
            
            % dependent variable - filtered the same way as above!
            this_data = all_data.s_all.(to_plot_2d{pp})(all_data.use_trial==1,dim_to_plot(dd));
            
            % now extract standard dev for each subj, condition
            
            % DV, IV, SUBJ
            thisX = nan(length(u_subj)*length(cu),3);
            cnt = 1;
            for ss = 1:length(u_subj)
                for cc = 1:length(cu)
                    thisidx = this_labels(:,1)==cu(cc) & this_labels(:,2)==ss;
                    thisX(cnt,:) = [std(this_data(thisidx)) cc ss];
                    cnt = cnt+1;
                    clear thisidx;
                end
            end
            
            allF(pp,dd,ii) = RMAOV1(thisX);
            
            for cp_idx = 1:size(cond_pairs,1)
%                 if ii == 1
%                     allT{cp_idx} = nan(length(to_plot_2d),length(dim_to_plot),niter+1);
%                 end
                [~,~,~,tmp_stats] = ttest(thisX(thisX(:,2)==cond_pairs(cp_idx,1),1),thisX(thisX(:,2)==cond_pairs(cp_idx,2),1));
                allT{cp_idx}(pp,dd,ii) = tmp_stats.tstat;
                clear tmp_stats;
            end
            
            clear cnt thisX;
            
        end
        
    end
end

fprintf('\n\n1-way repeated-measures ANOVA (against %i shuffling iterations):\n',niter);
for pp = 1:length(to_plot_2d)
    for dd = 1:length(dim_to_plot)
        fprintf('%s, %s:\tF = %0.03f, p = %0.03f\n',to_plot_2d{pp},dim_str{dd},allF(pp,dd,1),mean(squeeze(allF(pp,dd,2:end))>=allF(pp,dd,1)));
    end
end

fprintf('\n\nPaired t-test for eacn condition pair (against %i shuffling iterations):\n',niter);
for cp_idx = 1:size(cond_pairs,1)
    fprintf('\n%s vs %s\n',cond_str{cond_pairs(cp_idx,1)},cond_str{cond_pairs(cp_idx,2)});
    for pp = 1:length(to_plot_2d)
        for dd = 1:length(dim_to_plot)
            thisp = 2*mean( squeeze(abs(allT{cp_idx}(pp,dd,2:end))) >= abs(allT{cp_idx}(pp,dd,1))   );
            fprintf('%s, %s:\tT = %0.03f, p = %0.03f\n',to_plot_2d{pp},dim_str{dd},allT{cp_idx}(pp,dd,1),thisp);
                
        end
    end
    %fprintf('\n');
end