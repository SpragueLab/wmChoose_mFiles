% wmPri_plotEyeData.mm
%
% general plotting script for wmPri - plots sacc metrics as a function
% of condition (high, low priority)
%
% TCS 8/17/2018

% 11: drift correction
% 12: calibration (note: for run-wise, we're not using this at all..
% 13: fixation break
% 20: no primary saccade detected
% 21: bad primary saccade (too small/short)
% 22: large error for primary saccade

root = wmPri_loadRoot;

subj = {'AY','CC','KD','MR','SF','SH','XL','JK','YK','MH','SM'};
sess = {{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2','wmPri3'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'}};
%subj = {'JK','MH','SM','YK'};
%subj = {'JK','SM','YK'}; %JK, SM, YK are anti-prioritizers, along w/ MR
%sess = {{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'}};

%WHICH_EXCL = [11 13 20 21 22]; % don't exclude trials w/ calibration failures for now...
WHICH_EXCL = [20 21 22]; % don't exclude trials w/ calibration failures for now...

%incl_subj = [1 1 1 0 1 1 1 1 1 1 1]; % whether or not to include a subj for averaging...
incl_subj = [1 1 1 1 1 1 1 1 1 1 1]; % whether or not to include a subj for averaging...
%incl_subj = [1 1 1];

subj_colors = [0 0 0; 0.7 0.7 0.7]; % colors for included (1) or non-included (2) subj


% for now, let's use cat_struct to load/concatenate all data...
all_subj = nan(1000*length(subj),1);
u_subj = unique(cellfun(@(s) s(1:2),subj,'uniformoutput',0));

TARG_ECC = 12;

niter = 1000;

ET_HZ = 500;
ET_MS = 1000/ET_HZ; % sampling rate for plotting timeseries

all_data = [];

startidx = 1;

for ss = 1:length(subj)
    for sessidx = 1:length(sess{ss})
%     fn = sprintf('%s/data/%s_wmPri_scored.mat',root,subj{ss});
%     fprintf('Loading trial information from %s\n',fn);
%     this_data = load(fn);
    
    fn = sprintf('%s/wmPri_behav/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx});
    fprintf('Loading scored eye data from %s\n',fn);
    this_scored = load(fn);
    
    this_data.s_all = this_scored.ii_sess;
    this_data.sess_all = sessidx;
    
    this_subj = ss;%find(strcmpi(u_subj,subj{ss}(1:2)));
    
    all_data = cat_struct(all_data,this_data);
    all_subj(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1))) = this_subj;
    
    startidx = startidx+size(this_scored.ii_sess.trialinfo,1);
    
    clear this_subj this_data;
    end
end

% let's try this pattern for now
all_subj = all_subj(1:(startidx-1));
all_data.subj_all = all_subj;

% determine which trials to include
% first, narrow based on saccade preprocessing/scoring exclusions
% (wmChoose_extractSaccadeData1.m)
all_data.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data.s_all.excl_trial, 'UniformOutput',false));

% drop trials with very short (< 100 ms) or very long RT (> 1 s)
all_data.use_trial(all_data.s_all.i_sacc_rt<0.1 | all_data.s_all.i_sacc_rt>1.0) = 0;
%all_data.use_trial(all_data.s_all.i_sacc_err>5) = 0; % kill 'bad' trials (errors)

%% first, plot mean i_sacc, f_sacc error as a function of condition

mean_fig = figure;
scatter_fig = figure;

%to_plot = {'i_sacc_err','f_sacc_err','i_sacc_rt'};
to_plot = {'f_sacc_err','i_sacc_rt'};

cu = unique(all_data.s_all.trialinfo(:,1));

cond_str = {'High','Low'};

tmp_colors = lines(7);
cond_colors = tmp_colors([1 4],:);%lines(length(cu));
cond_pairs = [1 2];
%cond_pairs = [1 2; 2 3; 1 3]; % x, y axes of scatterplot


for pp = 1:length(to_plot)
    figure(mean_fig);
    subplot(1,length(to_plot),pp); hold on;
    
    thisd = nan(length(u_subj),length(cu));
    for cc = 1:length(cu)
        for ss = 1:length(u_subj)
            thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1;
            thisd(ss,cc) = mean(all_data.s_all.(to_plot{pp})(thisidx));
        end
        plot(cc+[-0.35 0.35],[1 1]*mean(thisd(:,cc)),'-','LineWidth',2.5,'Color',cond_colors(cc,:))
    end
    plot(1:length(cu),thisd.','-','Color',[0.5 0.5 0.5]);
    
    set(gca,'XTick',1:length(cu),'TickDir','out','LineWidth',1.5,'XTickLabel',cond_str,'FontSize',14,'XTickLabelRotation',-45);
    xlim([0.5 0.5+length(cu)]);
    title(to_plot{pp},'Interpreter','none');
    
end
   

figure(scatter_fig);
for pp = 1:length(to_plot)
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


% do a figure w/ error bars too....
figure;

stats_params_pval = nan(length(to_plot),1);
stats_params_tval = nan(length(to_plot),1);


for pp = 1:length(to_plot)
    subplot(1,length(to_plot),pp); hold on;
    
    thisd = nan(length(u_subj),length(cu));
    for cc = 1:length(cu)
        for ss = 1:length(u_subj)
            thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1;
            thisd(ss,cc) = mean(all_data.s_all.(to_plot{pp})(thisidx));
        end
    end
    
    % do stats
    [~,stats_params_pval(pp),~,tmpstats] = ttest(thisd(incl_subj==1,1),thisd(incl_subj==1,2));
    stats_params_tval(pp) = tmpstats.tstat; clear tmpstats;
    
    thisde = thisd - mean(thisd,2); % for within-subj SEM
    
    for cc = 1:length(cu)
        
        %plot(cc+[-0.35 0.35],[1 1]*mean(thisd(:,cc)),'-','LineWidth',2.5,'Color',cond_colors(cc,:))
        %thise = std(thisd(incl_subj==1,:),[],1)/sqrt(sum(incl_subj));
        thise = std(thisde(incl_subj==1,cc),[],1)/sqrt(sum(incl_subj));
        thism = mean(thisd(incl_subj==1,cc),1);
        
        plot(cc*[1 1],[-1 1]*thise+thism,'-','LineWidth',1.5,'Color',cond_colors(cc,:));
        plot(cc ,thism,'o','LineWidth',1.5,'Color',cond_colors(cc,:),'MarkerFaceColor','w','MarkerSize',10);
        clear thise thism;
    end
    for ss = 1:length(subj)
        if incl_subj(ss)==1
            plot(1:length(cu),thisd(ss,:).','-','Color',subj_colors(1,:));
        else
            plot(1:length(cu),thisd(ss,:).','-','Color',subj_colors(2,:));
        end
    end
    
    set(gca,'XTick',1:length(cu),'TickDir','out','XTickLabel',cond_str,'FontSize',14,'XTickLabelRotation',-45);
    xlim([0.5 0.5+length(cu)]);
    title(to_plot{pp},'Interpreter','none');
    
end


%% plot RT distribution for each subj, condition
figure;
for ss = 1:length(u_subj)
    for cc = 1:length(cu)
        
        subplot(length(cu),length(u_subj),(cc-1)*length(u_subj)+ss); hold on;
        
        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1;

        histogram(all_data.s_all.i_sacc_rt(thisidx),15,'BinLimits',[0 1.5],'FaceColor',cond_colors(cc,:));
        
        if ss == 1
            ylabel(cond_str{cc});
        end
        
        if cc == 1
            title(u_subj{ss});
        end
        
        if cc == length(cu)
            xlabel('RT (s)');
        end
        
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
            
            thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1;
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



%% radial saccade trajectories for each subj, condition (plotted going up/down)

% we care about EPOCH 5!!!

plot_dir = [1 -1]; % high, low priority

% compute R...
all_data.s_all.R = cellfun(@(x,y) sqrt(x.^2+y.^2),all_data.s_all.X,all_data.s_all.Y,'UniformOutput',false);

figure;
for ss = 1:length(subj)
    subplot(length(subj),1,ss); hold on;
    
    for cc = 1:length(cu)
        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1;
        thiscolor = cond_colors(cc,:); thisepochs = 5;
        %cellfun(@(x,d) plot((1:length(x))*ET_HZ, plot_dir(cc)*x(ismember(d,thisepochs)),'LineWidth',1,'Color',cond_colors(cc,:)),{all_data.s_all.X{thisidx}},{all_data.s_all.XDAT{thisidx}});
        cellfun(@(x,d) plot((1:sum(ismember(d,thisepochs))).'*ET_MS,plot_dir(cc)*x(ismember(d,thisepochs)),'LineWidth',1,'Color',cond_colors(cc,:)),{all_data.s_all.R{thisidx}},{all_data.s_all.XDAT{thisidx}});
        
    end
    
    if ss == length(subj)
        set(gca,'XTick',[0:500:2000],'TickDir','out');
        xlabel('Time (ms)');
    else
        set(gca,'XTick',[0:500:2000],'TickDir','out','XTickLabel',[]);
    end
    
    ylim([-17 17]); xlim([0 1500]);
    
    ylabel(subj{ss});
    set(gca,'YTick',[-12 -6 0 6 12],'YTickLabel',{'Low',[],[],[],'High'});
    
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
all_mu  = nan(length(dim_to_plot),length(to_plot_2d),length(cu),length(u_subj)); % store mean for each condition; subj (relative to target)

fig_dist = figure;
for pp = 1:length(to_plot_2d)
    for dd = 1:length(dim_to_plot)
        
        thisrow = (pp-1)*length(dim_to_plot)+dd;

        %figs_dist(pp,dd) = figure;
        
        for ss = 1:length(u_subj)
            subplot(this_nrows,length(u_subj),ss+(thisrow-1)*length(u_subj)); hold on;

            for cc = 1:length(cu)
                
                
                %subplot(length(cu),length(u_subj),ss+(cc-1)*length(u_subj)); hold on;
                
                
                thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1;
                
                [thish,thise] = histcounts(all_data.s_all.(to_plot_2d{pp})(thisidx,dim_to_plot(dd)) - dim_targ(dd),'BinWidth',0.5,'Normalization','pdf');
                plot(mean([thise(1:end-1);thise(2:end)],1),thish,'-','LineWidth',1,'Color',cond_colors(cc,:));
                
                all_mu(dd,pp,cc,ss) = mean(all_data.s_all.(to_plot_2d{pp})(thisidx,dim_to_plot(dd)) - dim_targ(dd));
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



% TODO: add a 3rd dimension of all_std which is the error collapsed across
% radial/tangential
%
% then, automatically plot and do stats on that dimension below

fig_std_sum = nan(length(dim_to_plot),1);

% cell: dimension (radial, tangential)
% value: saccade parameter (i_sacc, f_sacc)
stats_dim_tval{1} = nan(length(to_plot_2d),1);
stats_dim_tval{2} = nan(length(to_plot_2d),1);

stats_dim_pval{1} = nan(length(to_plot_2d),1);
stats_dim_pval{2} = nan(length(to_plot_2d),1);


% now plot all the std devs
for dd = 1:length(dim_to_plot)
    fig_std_sum(dd) = figure;

    for pp = 1:length(to_plot_2d)
        
        %subplot(length(to_plot_2d),length(dim_to_plot),dd+(pp-1)*length(to_plot_2d)); hold on;
        
        subplot(1,length(to_plot_2d),pp); hold on;
        
        for ss = 1:length(subj)
            if incl_subj(ss)==1
                plot(1:length(cu),squeeze(all_std(dd,pp,:,ss)),'-','Color',subj_colors(1,:));
            else
                plot(1:length(cu),squeeze(all_std(dd,pp,:,ss)),'-','Color',subj_colors(2,:));
            end
        end
        
        
        thisd = squeeze(all_std(dd,pp,:,:)).'; % subj x cond
        thisde = thisd - mean(thisd,2);
        
        [~,stats_dim_pval{dd}(pp),~,tmpstats] = ttest(thisd(incl_subj==1,1),thisd(incl_subj==1,2));
        stats_dim_tval{dd}(pp) = tmpstats.tstat; clear tmpstats;
        
        
        for cc = 1:length(cu)
            
            thismu = mean(squeeze(all_std(dd,pp,cc,incl_subj==1)));
            %thise  = std( squeeze(all_std(dd,pp,cc,incl_subj==1)))/sqrt(sum(incl_subj));
            
            % within-subj SEM
            thise  = std( squeeze(thisde(incl_subj==1,cc)))/sqrt(sum(incl_subj));
            
            
            %plot(cc+[-0.35 0.35],thismu*[1 1],'-','LineWidth',2,'Color',cond_colors(cc,:));
            plot(cc*[1 1],thismu+[-1 1]*thise,'-','LineWidth',1.5,'Color',cond_colors(cc,:));
            plot(cc,thismu,'o','LineWidth',1.5,'MarkerSize',10,'MarkerFaceColor','w','Color',cond_colors(cc,:));
            
            
            clear thismu;
        end
        
        set(gca,'XTick',1:length(cu),'XTickLabel',cond_str,'XTickLabelRotation',-45,'TickDir','out','FontSize',14);
        xlim([0.5 length(cu)+0.5]);
        
        %if dd == 1
            title(to_plot_2d{pp},'Interpreter','none');
        %end
        
        if pp == 1
            ylabel('Error (std dev; \circ)');
        end
        
    end
    match_ylim(get(fig_std_sum(dd),'Children'));
    set(gcf,'Name',sprintf('%s - standard deviation',dim_str{dd}),'NumberTitle','off');

    %sgtitle(dim_str{dd});
    
end


%% also do this for 'average' error (avg of rad/tang std dev)

fig_std_sum_avg = figure;
for pp = 1:length(to_plot_2d)
    
        
        subplot(length(to_plot_2d),1,pp); hold on;
        
        plot(1:length(cu),squeeze(mean(all_std(:,pp,:,:),1)),'o-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w');
        
        for cc = 1:length(cu)
            
            thismu = mean(squeeze(mean(all_std(:,pp,cc,:),1)));
            thise  = std( squeeze(mean(all_std(:,pp,cc,:),1)),[],1)/sqrt(length(subj));
            %plot(cc+[-0.35 0.35],thismu*[1 1],'-','LineWidth',2,'Color',cond_colors(cc,:));
            plot(cc*[1 1],[-1 1]*thise+thismu,'-','LineWidth',1.5,'Color',cond_colors(cc,:));
            plot(cc,thismu,'o','LineWidth',1.5,'Color',cond_colors(cc,:),'MarkerSize',7,'MarkerFaceColor','w');
            
            clear thismu;
        end
        
        set(gca,'XTick',1:length(cu),'XTickLabel',cond_str,'XTickLabelRotation',-45);
        xlim([0.35 length(cu)+0.65]);
        
        
        ylabel(to_plot_2d{pp},'Interpreter','none');
        
        
        if pp == 1
            title('Average error (both dimensions)');
        end
        
        set(gca,'XTick',1:length(cu),'TickDir','out','LineWidth',1.5,'XTickLabel',cond_str,'FontSize',14,'XTickLabelRotation',-45);
        xlim([0.5 0.5+length(cu)]);
    
end
match_ylim(get(fig_std_sum_avg,'Children'));
set(gcf,'Name','Standard deviation - average','NumberTitle','off');





%% and all the means (bias)
fig_mu_sum = figure;
for pp = 1:length(to_plot_2d)
    for dd = 1:length(dim_to_plot)
        
        subplot(length(to_plot_2d),length(dim_to_plot),dd+(pp-1)*length(to_plot_2d)); hold on;
        
        plot(1:length(cu),squeeze(all_mu(dd,pp,:,:)),'o-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w');
        
        for cc = 1:length(cu)
            
            thismu = mean(squeeze(all_mu(dd,pp,cc,:)));
            
            plot(cc+[-0.35 0.35],thismu*[1 1],'-','LineWidth',2,'Color',cond_colors(cc,:));
            
            clear thismu;
        end
        
        set(gca,'XTick',1:length(cu),'XTickLabel',cond_str,'XTickLabelRotation',-45);
        xlim([0.35 length(cu)+0.65]);
        
        if dd == 1
            ylabel(to_plot_2d{pp},'Interpreter','none');
        end
        
        if pp == 1
            title(dim_str{dd});
        end
        
    end
end
match_ylim(get(fig_mu_sum,'Children'));
set(gcf,'Name','Bias','NumberTitle','off');

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


% ~~~~~ for now, hack and just do dd==3 on its own...

allF = nan(length(to_plot_2d),niter+1);   % params x dimensions x iterations?
allT = cell(size(cond_pairs,1),1);  % each of these same dims as above
for cp_idx = 1:size(cond_pairs,1)
    allT{cp_idx} = nan(length(to_plot_2d),niter+1);
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
        %for dd = 1:length(dim_to_plot)
            
            % dependent variable - filtered the same way as above!
            this_data = all_data.s_all.(to_plot_2d{pp})(all_data.use_trial==1,:);
            
            % now extract standard dev for each subj, condition
            
            % DV, IV, SUBJ
            thisX = nan(length(u_subj)*length(cu),3);
            cnt = 1;
            for ss = 1:length(u_subj)
                for cc = 1:length(cu)
                    thisidx = this_labels(:,1)==cu(cc) & this_labels(:,2)==ss;
                    thisX(cnt,:) = [mean(std(this_data(thisidx),[],1),2) cc ss];
                    cnt = cnt+1;
                    clear thisidx;
                end
            end
            
            allF(pp,ii) = RMAOV1(thisX);
            
            for cp_idx = 1:size(cond_pairs,1)
%                 if ii == 1
%                     allT{cp_idx} = nan(length(to_plot_2d),length(dim_to_plot),niter+1);
%                 end
                [~,~,~,tmp_stats] = ttest(thisX(thisX(:,2)==cond_pairs(cp_idx,1),1),thisX(thisX(:,2)==cond_pairs(cp_idx,2),1));
                allT{cp_idx}(pp,ii) = tmp_stats.tstat;
                clear tmp_stats;
            end
            
            clear cnt thisX;
            
        %end
        
    end
end

fprintf('\n\nAVERAGE ERROR (average of radial/tangential std dev\n');
fprintf('1-way repeated-measures ANOVA (against %i shuffling iterations):\n',niter);
for pp = 1:length(to_plot_2d)
    %for dd = 1:length(dim_to_plot)
        fprintf('%s, avg:\tF = %0.03f, p = %0.03f\n',to_plot_2d{pp},allF(pp,1),mean(squeeze(allF(pp,2:end))>=allF(pp,1)));
    %end
end

fprintf('\n\nPaired t-test for eacn condition pair (against %i shuffling iterations):\n',niter);
for cp_idx = 1:size(cond_pairs,1)
    fprintf('\n%s vs %s\n',cond_str{cond_pairs(cp_idx,1)},cond_str{cond_pairs(cp_idx,2)});
    for pp = 1:length(to_plot_2d)
        %for dd = 1:length(dim_to_plot)
            thisp = 2*mean( squeeze(abs(allT{cp_idx}(pp,2:end))) >= abs(allT{cp_idx}(pp,1))   );
            fprintf('%s, avg:\tT = %0.03f, p = %0.03f\n',to_plot_2d{pp},allT{cp_idx}(pp,1),thisp);
                
        %end
    end
    %fprintf('\n');
end
