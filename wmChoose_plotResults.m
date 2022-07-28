% wmChoose_plotResults.m
%
% general plotting script for wmChoose - plots sacc metrics as a function
% of condition (R1, R2-cued, R2-choose)
%
% TCS 4/13/2018

% 11: drift correction
% 12: calibration (note: for run-wise, we're not using this at all..
% 13: fixation break
% 20: no primary saccade detected
% 21: bad primary saccade (too small/short)
% 22: large error for primary saccade

root = 'Z:/projects/wmChoose';
subj = {'sub001','sub002','sub003','sub004','sub005','sub006','sub007','sub008','sub009','sub010','sub011','sub012','sub013','sub014','sub016','sub017','sub018','sub019','sub020','sub021'};
%subj = {'sub001','sub003','sub004','sub005','sub007','sub008','sub010','sub011','sub012','sub013','sub016','sub017','sub018','sub019','sub020'};

%WHICH_EXCL = [11 13 20 21 22]; % don't exclude trials w/ calibration failures for now...
WHICH_EXCL = [13 20 21]; % don't exclude trials w/ calibration failures for now...

% for now, let's use cat_struct to load/concatenate all data...
all_subj = nan(1000*length(subj),1);
%u_subj = unique(cellfun(@(s) s(1:2),subj,'uniformoutput',0));
u_subj = subj;

TARG_ECC = 12;

niter = 100; % we don't use this much now - if we do, set to 1000


ET_HZ = 1000;
ET_MS = 1000/ET_HZ; % sampling rate for plotting timeseries

all_data = [];

startidx = 1;

for ss = 1:length(subj)
    fn = sprintf('%s/data/%s_wmChoose_behav.mat',root,subj{ss});
    fprintf('Loading trial information from %s\n',fn);
    this_data = load(fn);
    
    fn = sprintf('%s/data/%s_wmChoose_scored.mat',root,subj{ss});
    fprintf('Loading scored eye data from %s\n',fn);
    this_scored = load(fn);
    
    this_data.s_all = this_scored.ii_sess;
    
    this_subj = ss; %find(strcmpi(u_subj,subj{ss}(1:2)));
    
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


%% first, plot i_sacc, f_sacc error as a function of condition

mean_fig = figure;
scatter_fig = figure;

%to_plot = {'i_sacc_err','f_sacc_err','i_sacc_rt'};
to_plot = {'f_sacc_err','i_sacc_rt'};
cu = unique(all_data.c_all(:,1));

%cond_str = {'R1','R2-cue','R2-choose'};
cond_str = {'R1','R2-random','R2-best'};

cond_colors = lines(length(cu));

cond_pairs = [1 2; 2 3; 1 3]; % x, y axes of scatterplot


for pp = 1:length(to_plot)

    fprintf('\nParametric stats for %s\n',to_plot{pp});

    figure(mean_fig);
    subplot(1,length(to_plot),pp); hold on;
    
    thisd = nan(length(u_subj),length(cu));
    thisd_anova = nan(length(u_subj)*length(cu),1);
    thisX_anova = nan(length(u_subj)*length(cu),2); % cond,subj
    cnt = 1;
    for cc = 1:length(cu)
        for ss = 1:length(u_subj)
            thisidx = all_data.subj_all==ss & all_data.c_all(:,1)==cu(cc) & all_data.use_trial==1;
            thisd(ss,cc) = mean(all_data.s_all.(to_plot{pp})(thisidx));  %change of std or mean here
            thisd_anova(cnt) = thisd(ss,cc);
            thisX_anova(cnt,:) = [cu(cc) ss];
            cnt = cnt+1;
        end
        %plot(cc+[-0.35 0.35],[1 1]*mean(thisd(:,cc)),'-','LineWidth',2.5,'Color',cond_colors(cc,:))
    end
    plot(1:length(cu),thisd.','-','Color',[0.5 0.5 0.5]);
    for cc = 1:length(cu)
        %for ss = 1:length(u_subj)
        plot(cc,thisd(:,cc),'o','MarkerSize',1.5,'Color',cond_colors(cc,:),'MarkerFaceColor',cond_colors(cc,:));

        thism = mean(thisd(:,cc)); thise = std(thisd(:,cc))/sqrt(length(subj));

        plot(cc*[1;1],thism+[-1;1]*thise,'-','Color',cond_colors(cc,:),'LineWidth',2);
        plot(cc,mean(thisd(:,cc)),'o','MarkerSize',7,'LineWidth',2,'Color',cond_colors(cc,:),'MarkerFaceColor','w');

        %end
    end
    %errorbar(1:length(cu),mean(thisd),std(thisd)/sqrt(length(thisd)),'--ko','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black','CapSize',10);

    set(gca,'XTick',1:length(cu),'TickDir','out','LineWidth',1,'XTickLabel',cond_str,'FontSize',14,'XTickLabelRotation',-45);
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

        [~,thisp,~,thisstats] = ttest(thisd(:,cond_pairs(cp,1)),thisd(:,cond_pairs(cp,2)));

        fprintf('T-test: %s vs %s, T(%i) = %.03f, p = %0.05f, dz = %0.05f\n',cond_str{cond_pairs(cp,1)},cond_str{cond_pairs(cp,2)}, thisstats.df,thisstats.tstat,thisp,thisstats.tstat/sqrt(length(subj)));

        clear thisp thisstats;

    end


    % stats for this parameter (ANOVA)
    [thisp,thisanovatab] = anovan( thisd_anova,thisX_anova  , 'random',2,'varnames',{'condition','subj'},'model','interaction' ,'display','off');

    % compute partial eta squared (SS_cond/(SS_cond+SS_err))
    this_peta2 = thisanovatab{2,2}/(thisanovatab{2,2}+thisanovatab{4,2}); % SS_cond/(SS_cond + SS_cond*subj)
    
    %[thisp,thisanovatab] = anovan( thisd_anova,thisX_anova  , 'random',2,'varnames',{'condition','subj'},'model','interaction');
    fprintf('1-way RM AOV: F(%i,%i) = %.03f, p = %0.05f, partial eta^2 = %0.05f\n\n',thisanovatab{2,3},thisanovatab{4,3},thisanovatab{2,6},thisp(1),this_peta2);

    clear thisp thisanovatab thisd_anova thisX_anova;
    
end



%% plot RT distribution for each subj, condition
figure;
for ss = 1:length(u_subj)
    for cc = 1:length(cu)
        
        subplot(length(cu),length(u_subj),(cc-1)*length(u_subj)+ss); hold on;
        
        thisidx = all_data.subj_all==ss & all_data.c_all(:,1)==cu(cc) & all_data.use_trial==1;

        histogram(all_data.s_all.i_sacc_rt(thisidx),45,'BinLimits',[0 1.5],'FaceColor',cond_colors(cc,:));
        
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



%% for grants, plot just R2-cued and R2-choose, and just mean +- SEM
if 0 % typically don't need to plot this figure
figure; hold on;
cu_grant = [2 3];
tmpd = nan(length(subj),length(cu_grant));

for cc = 1:length(cu_grant)
    for ss = 1:length(subj)
        thisidx = all_data.subj_all==ss & all_data.c_all(:,1)==cu_grant(cc) & all_data.use_trial==1;
        tmpd(ss,cc) = mean(all_data.s_all.('f_sacc_err')(thisidx));
    end
    thise = std(tmpd(:,cc),[],1)./sqrt(length(subj));
    thism = mean(tmpd(:,cc),1);

    plot(cc*[1 1],thism+[-1 1].*thise,'-','LineWidth',1,'Color',cond_colors(cu_grant(cc),:));
    plot(cc,thism,'o','MarkerFaceColor','w','MarkerSize',10','LineWidth',1,'Color',cond_colors(cu_grant(cc),:));
end

ylim([1.3 2.3]);xlim([.5 2.5]);
set(gca,'LineWidth',1,'TickDir','out','XTick',[1 2],'XTickLabel',{'Cued','Chosen'},'YTick',[1.25:.25:2.5],'FontSize',14);
ylabel('Error');
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
                
                
                thisidx = all_data.subj_all==ss & all_data.c_all(:,1)==cu(cc) & all_data.use_trial==1;
                
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
set(gcf,'Name','Standard deviation','NumberTitle','off');


%% also do this for 'average' error (avg of rad/tang std dev)

fig_std_sum_avg = figure;
for pp = 1:length(to_plot_2d)
    
        
        subplot(length(to_plot_2d),1,pp); hold on;
        
        plot(1:length(cu),squeeze(mean(all_std(:,pp,:,:),1)),'o-','Color',[0.5 0.5 0.5]);
        
        for cc = 1:length(cu)
            
            thismu = mean(squeeze(mean(all_std(:,pp,cc,:),1)));
            
            plot(cc+[-0.35 0.35],thismu*[1 1],'-','LineWidth',2,'Color',cond_colors(cc,:));
            
            clear thismu;
        end
        
        set(gca,'XTick',1:length(cu),'XTickLabel',cond_str,'XTickLabelRotation',-45);
        xlim([0.35 3.65]);
        
        if dd == 1
            ylabel(to_plot_2d{pp},'Interpreter','none');
        end
        
        if pp == 1
            title('Average error (both dimensions)');
        end
        
    
end
match_ylim(get(fig_std_sum_avg,'Children'));
set(gcf,'Name','Standard deviation - average','NumberTitle','off');



% and all the means
fig_mu_sum = figure;
for pp = 1:length(to_plot_2d)
    for dd = 1:length(dim_to_plot)
        
        subplot(length(to_plot_2d),length(dim_to_plot),dd+(pp-1)*length(to_plot_2d)); hold on;
        
        plot(1:length(cu),squeeze(all_mu(dd,pp,:,:)),'o-','Color',[0.5 0.5 0.5]);
        
        for cc = 1:length(cu)
            
            thismu = mean(squeeze(all_mu(dd,pp,cc,:)));
            
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
match_ylim(get(fig_mu_sum,'Children'));
set(gcf,'Name','Bias','NumberTitle','off');

% TODO: 
% - trial exclusion % by condition (maybe broken down by exclusion type?)
%% table of % included trials from each subj

percent_useTrial = nan(length(u_subj),1);
percent_excTrial = nan(length(u_subj),1);
for cc = 1:length(cu)
    for ss = 1:length(u_subj)

        thisidx = all_data.subj_all==ss & all_data.c_all(:,1)==cu(cc) & all_data.use_trial==1;
        thisidx2 = all_data.subj_all==ss & all_data.c_all(:,1)==cu(cc) & all_data.use_trial==0;
        thisidx3 = all_data.subj_all==ss & all_data.c_all(:,1)==cu(cc);
        percent_useTrial(ss) = sum(thisidx)/sum(thisidx3);
        percent_excTrial(ss) = sum(thisidx2)/sum(thisidx3);
    end

end

% - trace for each condition; distribution for each condition

%% response saccade trajectories for examplar subj, condition

% we care about EPOCH 3 here!!!

plot_dir = [1 1 1]; % R1, R2v, R2c
plot_v_offset = [30 15 0];

% compute R...
all_data.s_all.R = cellfun(@(x,y) sqrt(x.^2+y.^2),all_data.s_all.X,all_data.s_all.Y,'UniformOutput',false);

which_subj_traj = 1; % to plot all, do which_subj_traj = subj;

figure;
for ss = 1:length(which_subj_traj)
    subplot(length(which_subj_traj),1,ss); hold on;
    
    for cc = 1:length(cu)
        thisidx = all_data.subj_all==which_subj_traj(ss) & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1;
        thiscolor = cond_colors(cc,:); thisepochs = [3];
        %cellfun(@(x,d) plot((1:length(x))*ET_HZ, plot_dir(cc)*x(ismember(d,thisepochs)),'LineWidth',1,'Color',cond_colors(cc,:)),{all_data.s_all.X{thisidx}},{all_data.s_all.XDAT{thisidx}});
        cellfun(@(x,d) plot((1:sum(ismember(d,thisepochs))).'*ET_MS,plot_v_offset(cc)+plot_dir(cc)*x(ismember(d,thisepochs)),'LineWidth',1,'Color',cond_colors(cc,:)),{all_data.s_all.R{thisidx}},{all_data.s_all.XDAT{thisidx}});
        
    end
    
    if ss == length(which_subj_traj)
        set(gca,'XTick',[0:500:1500],'TickDir','out');
        xlabel('Time (ms)');
    else
        set(gca,'XTick',[0:500:1500],'TickDir','out','XTickLabel',[]);
    end
    
    ylim([-3 45]); xlim([0 1500]);
    
    ylabel(subj{which_subj_traj(ss)});
    %set(gca,'YTick',[-12 -6 0 6 12],'YTickLabel',{'Low',[],[],[],'High'});
    set(gca,'YTick',[0 6 12 15 21 27 30 36 42],'YTickLabel',{[],'R2-choose',[],[],'R2-cue',[],[],'R1',[]});
    % TODO: do this programmatically....
    
end

%% Trajectory for each trial in one subject
figure;
for ss = 1:length(which_subj_traj)
    subplot(length(which_subj_traj),1,ss); hold on;
    
    for cc = 1:length(cu)
        thisidx = all_data.subj_all==which_subj_traj(ss) & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1 & all_data.r_all==5;
        thiscolor = cond_colors(cc,:); thisepochs = [3];
        %cellfun(@(x,d) plot((1:length(x))*ET_HZ, plot_dir(cc)*x(ismember(d,thisepochs)),'LineWidth',1,'Color',cond_colors(cc,:)),{all_data.s_all.X{thisidx}},{all_data.s_all.XDAT{thisidx}});
        %cellfun(@(x,d) plot((1:sum(ismember(d,thisepochs))).'*ET_MS,plot_v_offset(cc)+plot_dir(cc)*x(ismember(d,thisepochs)),'LineWidth',1,'Color',cond_colors(cc,:)),{all_data.s_all.R{thisidx}},{all_data.s_all.XDAT{thisidx}});
        thisidx2 = find(thisidx);

        for ii = 1:length(thisidx2)
            plot(all_data.s_all.X{thisidx2(ii)},all_data.s_all.Y{thisidx2(ii)},'Color',cond_colors(cc,:));
            plot(all_data.s_all.targ(thisidx2,1),all_data.s_all.targ(thisidx2,2),'ro','MarkerFaceColor','r','MarkerSize',6);
        end
    end

    plot(0,0,'k+','MarkerSize',8,'LineWidth',1.5);
   
    xlim([-15 15]);ylim([-15 15]);
    axis square

    tmpx = 12*cos(linspace(0,2*pi,1001));
    tmpy = 12*sin(linspace(0,2*pi,1001)); 
    plot(tmpx,tmpy,'k-','LineWidth',0.5);

    set(gca,'TickDir','out','XTick',[-12 0 12],'YTick',[-12 0 12]); 
%     if ss == length(which_subj_traj)
%         set(gca,'XTick',[0:500:1500],'TickDir','out');
%         xlabel('Time (ms)');
%     else
%         set(gca,'XTick',[0:500:1500],'TickDir','out','XTickLabel',[]);
%     end
%     
%     ylim([-3 45]); xlim([0 1500]);
%     
%     ylabel(subj{which_subj_traj(ss)});
%     %set(gca,'YTick',[-12 -6 0 6 12],'YTickLabel',{'Low',[],[],[],'High'});
%     set(gca,'YTick',[0 6 12 15 21 27 30 36 42],'YTickLabel',{[],'R2-choose',[],[],'R2-valid',[],[],'R1',[]});
%     % TODO: do this programmatically....
    
end


%% Target loc against Report loc from examplar subj
    
%compute target theta
all_data.s_all.TargetX = num2cell(all_data.s_all.targ(:,1));
all_data.s_all.TargetY = num2cell(all_data.s_all.targ(:,2));
%all_data.s_all.TarT = cellfun(@(x,y) atan2d(y,x), all_data.s_all.TargetX,all_data.s_all.TargetY,'UniformOutput', false);
all_data.s_all.TarT = cellfun(@(x,y) atan2d(y,x), all_data.s_all.TargetX,all_data.s_all.TargetY);
%all_data.s_all.TarT360 = cellfun(@(x) wrapTo360(x), all_data.s_all.TarT,'UniformOutput', false);
%all_data.s_all.TarT360 = cellfun(@(x) wrapTo360(x), all_data.s_all.TarT);
all_data.s_all.TarT360 = mod(all_data.s_all.TarT+360,360);


%compute reported theta
all_data.s_all.f_saccX = num2cell(all_data.s_all.f_sacc_raw(:,1));
all_data.s_all.f_saccY = num2cell(all_data.s_all.f_sacc_raw(:,2));
all_data.s_all.RepT = cellfun(@(x,y) atan2d(y,x), all_data.s_all.f_saccX,all_data.s_all.f_saccY);
all_data.s_all.RepT360 = mod(all_data.s_all.RepT+360,360);

figure;
for ss = 1:length(u_subj)

    
    
    for cc = 1:length(cu)

        subplot(length(cu),length(u_subj),ss+length(u_subj)*(cc-1)); hold on;

        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1;
        thisidx2 = find(thisidx);
        

        %for ii = 1:length(thisidx2)    
           
           scatter(all_data.s_all.TarT360(thisidx2),all_data.s_all.RepT360(thisidx2),10,'r','filled', 'markerfacealpha',.5);
           
        %end

     xlim([0 360]);ylim([0 360]);
     set(gca,'TickDir','out','XTick',[0 180 360],'YTick',[0 180 360]); 
     xlabel('Target location (°)'); ylabel('Reported location (°)');
     axis square;
    end



end

figure;
for ss = 1:length(which_subj_traj)

    
    
    for cc = 1:length(cu)

        %subplot(length(cu),length(u_subj),ss+length(u_subj)*(cc-1)); hold on;

        thisidx = all_data.subj_all==which_subj_traj(ss) & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1;
        thisidx2 = find(thisidx);
        

        %for ii = 1:length(thisidx2)    
           
           scatter(all_data.s_all.TarT360(thisidx2),all_data.s_all.RepT360(thisidx2),25,'r','filled', 'markerfacealpha',.5);
           
        %end

     xlim([0 360]);ylim([0 360]);
     set(gca,'TickDir','out','XTick',[0 180 360],'YTick',[0 180 360]); 
     xlabel('Target location (°)'); ylabel('Reported location (°)');
     axis square;
    end

end

%% plot error distribution for examplar subj

all_data.s_all.f_sacc_180pol = all_data.s_all.TarT - all_data.s_all.RepT;

figure;
 
for ss = 1:length(which_subj_traj)
    
    
    for cc = 1:length(cu)
        thisidx = all_data.subj_all==which_subj_traj(ss) & all_data.s_all.trialinfo(:,1)==cu(1) & all_data.use_trial==1;
        
        thisidx2 = find(thisidx);
  
          
            histogram(all_data.s_all.f_sacc_180pol(thisidx2),10);
  
        
    end

    xlim([-180 180]);ylim([0 100]);
    set(gca,'TickDir','out','XTick',[-180 -90 0 90 180],'YTick',[0 50 100]); 
    xlabel('Behavioral error (°)'); ylabel('Trial count');

end



%% plot polarhistogram of target location of choose condition for each subj

%compute radian
all_data.s_all.TarRadian = deg2rad(all_data.s_all.TarT360);

figure;
for ss = 1:length(u_subj)
    for cc = 1:length(cu)

        subplot(length(cu),length(u_subj),(cc-1)*length(u_subj)+ss); hold on;

        thisidx = all_data.subj_all==ss & all_data.c_all(:,1)==cu(cc) & all_data.use_trial==1;
   
           
           histogram(all_data.s_all.TarRadian(thisidx),16);

           if cc == 1
            title(sprintf('%s-%s',u_subj{ss},cond_str{cc}));
           end
        

    end

end


%% seperate split half session for target loc choice for choose condition

figure;
for ss = 1:length(u_subj)

    %his_split_half = [];
    tmpr = all_data.r_all(all_data.subj_all==ss & all_data.c_all(:,1)==3);
    nr = length(unique(tmpr)); % number of runs for this subj

    this_split_half = (tmpr>ceil(nr/2))+1;

    for ii = 1:2

        subplot(2,length(u_subj),(ii-1)*length(u_subj)+ss); hold on;

        thisidx = all_data.subj_all==ss & all_data.c_all(:,1)==3;% & all_data.use_trial==1;
        tmpd = all_data.s_all.TarRadian(thisidx);
        
        tmpuse = all_data.use_trial(thisidx);
        
        tmpd = tmpd(this_split_half==ii & tmpuse == 1);

           histogram(tmpd,8);
 
            if ss ~=1
                set(gca,'YTickLabel',[]);
            end

           if cc == 1
            title(sprintf('%s-%s',u_subj{ss},cond_str{cc}));
           end
        

    end

end

match_ylim(get(gcf,'Children'));
%same_hemi_trials

%% target color investigation


% we want:
% - n_trials x 2 indices into set of unique colors
% - col 1: reported color index
% - col 2: non-reported color (or nan)

% first - figure out the unique colors

% 200 0 0;        % red.78 0 0
% 0 0 255;        % blue 0 0 1
% 180 0 180;      % purple .71 0 .71
% 130 130 0;      % yellow .51 .51 0

u_colors = unique(all_data.colors_all{20,1},'rows');% get the unique value from the data set see how many color stimuli we had
u_colors_rgb = [0 0 1;
    .51 .51 0;
    .71 0 .71;
    .78 0 0];
%u_color_rgb = rgb2hsv(u_colors_rgb);
all_targ_colors_idx = nan(length(all_subj),2);% create a data structure for color index, each sub from each trial will have two values for the two stimuli color, one for R1

% create the color index
for ss = 1:length(u_subj)

    this_info = all_data.s_all.trialinfo(all_subj == ss,:);
    this_colors_idx = nan(size(this_info,1),2);


    for tt = 1:size(this_info,1)
        this_targ = this_info(tt,6);
        if this_targ == 1
            this_colors_idx(tt,1) = find(ismember(u_colors,all_data.colors_all{ss,1}(tt,:),"rows"));

            if this_info(tt,1) ~= 1
                this_colors_idx(tt,2) = find(ismember(u_colors,all_data.colors_all{ss,2}(tt,:),"rows"));


            end
        else
            this_colors_idx(tt,1) = find(ismember(u_colors,all_data.colors_all{ss,2}(tt,:),"rows"));

            if this_info(tt,1) ~= 1
                this_colors_idx(tt,2) = find(ismember(u_colors,all_data.colors_all{ss,1}(tt,:),"rows"));


            end
        end
    end
    all_targ_colors_idx(all_subj == ss,:) = this_colors_idx;
end

%% plot accuracy as function of target color 
my_Colors = 1:1:4;

figure;
for ss = 1:length(u_subj)
    subplot(4,5,ss); hold on
    for tt = 1:length(my_Colors)

        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(1) & all_data.use_trial==1 & all_targ_colors_idx == tt;

        plot(my_Colors(tt),mean(all_data.s_all.f_sacc_err(thisidx)),'o','Color',u_colors_rgb(tt,:),'MarkerFaceColor',u_colors_rgb(tt,:));
    end
    xlim([0 5]);ylim([0 5]);
    xlabel('Target color'); ylabel('Final Saccade Error');
end

% same as above but collapses across subs

figure;
for ss = 1:length(u_subj)

    for tt = 1:length(my_Colors)
        subplot(1,1,1); hold on
        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(1) & all_data.use_trial==1 & all_targ_colors_idx == tt;

        plot(my_Colors(tt),mean(all_data.s_all.f_sacc_err(thisidx)),'o','Color',u_colors_rgb(tt,:),'MarkerFaceColor',u_colors_rgb(tt,:));
    end
    xlim([0 5]);ylim([0 5]);
    xlabel('Target color'); ylabel('Final Saccade Error');
end
%% plot R1 error against target color choice in R2Choose

all_col_corr = nan(length(u_subj),1);

figure;
for ss = 1:length(u_subj)% loop through each subj

     this_n_choose = sum(all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(3) & all_data.use_trial==1);% compute all number of trials in choose condition
    tmpErr = nan(length(my_Colors),1);
    tmpP = nan(length(my_Colors),1);


    for tt = 1:length(my_Colors)

        subplot(4,5,ss);hold on

        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(1) & all_data.use_trial==1 & all_targ_colors_idx(:,1) == tt;% index trials used in one object condition
        thisidx2 = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(3) & all_data.use_trial==1 & all_targ_colors_idx(:,1) == tt;% index trials used in choose condition
        this_x = mean(all_data.s_all.f_sacc_err(thisidx));
        this_y = sum(thisidx2)/this_n_choose;
        plot(this_x,this_y,'o','Color',u_colors_rgb(tt,:),'MarkerFaceColor',u_colors_rgb(tt,:));


        tmpErr(tt) = this_x;
        tmpP (tt) = this_y;
    end

    all_col_corr(ss) = corr(tmpErr,tmpP);

    xlim([0 5]);ylim([0 .5]);
    set(gca,'TickDir','out','XTick',[0 1 5],'YTick',[0 0.1 .5]);
    xlabel('R1 Error'); ylabel('Proportion of Choose Color in R2 Choose');

end
match_xlim(get(gcf,'Children')); match_ylim(get(gcf,'Children'));

figure; hold on;

plot(ones(size(all_col_corr)),all_col_corr,'ko');
plot(1,mean(all_col_corr),'ko','MarkerSize',15,'LineWidth',1.5,'MarkerFaceColor','k');

all_col_corr_z = atanh(all_col_corr);

% compute and spit out stats for this test
[~,thisp,~,thisstats] = ttest(all_col_corr_z);

fprintf('\nColor heuristic analysis: t-test of correlations against 0\n');
fprintf('T-test: corr against zero, T(%i) = %.05f, p = %.05f, dz = %.05f\n\n',thisstats.df,thisstats.tstat,thisp,thisstats.tstat/sqrt(length(subj)));


clear thisp thisstats;

%% plot accuracy as function of location bin

%comment each line of the below two graphs to make sure I can fully
%understand

all_data.s_all.f_sacc_y = all_data.s_all.f_sacc(:,2);

my_Bin = 0:30:360; % create bin limit for 0 to 360, make 1 bin for 30 degree interval therefore result in 12 bins in total
%[BINS, EDGES] = discretize(X,N) where N is a scalar integer, divides the range of X into N uniform bins, and also returns the bin edges.
location_Bin = discretize(all_data.s_all.TarT360,my_Bin);%put target location data into bin from 1-12, e.g., 1 is 0-30 degree and so on.

bin_colors = hsv(length(my_Bin));%give each bin a distinct color

figure;
for ss = 1:length(u_subj) % go through all the subjects


    subplot(4,5,ss);hold on % create 20 plots in 4 by 5 arrangement
    for bb = 1:(length(my_Bin)-1) % loop over 1 to 11 because we want 12 interval from 13 numbers, 0 and 360 is the same here

        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(1) & all_data.use_trial==1 & location_Bin==bb;
        thisidx2 = find(thisidx);% minimize the data size for plot by selecting a group of data with the above criterions


        plot((my_Bin(bb+1)+my_Bin(bb))/2,mean(all_data.s_all.f_sacc_err(thisidx2)),'o','Color',bin_colors(bb,:),'MarkerFaceColor',bin_colors(bb,:),'MarkerSize',4);
    end %plot(data point on x axis is ploted on the center of each interval e.g., (0+30)/2 = 15 and so on for all bin

    text(330,4.5,u_subj{ss},'FontAngle','italic','HorizontalAlignment','right');

    xlim([0 360]);ylim([0 5]);
    xticks([0 90 180 270 360]);
    yticks(0:2:6);
    % set(gca,'TickDir','out','XTick',[-180 0 180],'YTick',[-180 0 180]);
    if ss == 16; xlabel('Target location (°)'); ylabel('Final Saccade Error (MAE; °)'); end
    set(gca,'TickDir','out','LineWidth',0.75);

end

%same graph collapse across sub
figure;
for ss = 1:length(u_subj) % go through all the subjects
    for bb = 1:(length(my_Bin)-1) % loop over 1 to 11 because we want 12 interval from 13 numbers, 0 and 360 is the same here
        subplot(1,1,1);hold on
        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(1) & all_data.use_trial==1 & location_Bin==bb;
        thisidx2 = find(thisidx);% minimize the data size for plot by selecting a group of data with the above criterions

        plot((my_Bin(bb+1)+my_Bin(bb))/2,mean(all_data.s_all.f_sacc_err(thisidx2)),'o','Color',bin_colors(bb,:),'MarkerFaceColor',bin_colors(bb,:));
        
    end %plot(data point on x axis is ploted on the center of each interval e.g., (0+30)/2 = 15 and so on for all bin

    xlim([0 360]);ylim([0 5]);
    xticks([0 90 180 270 360]);
    % set(gca,'TickDir','out','XTick',[-180 0 180],'YTick',[-180 0 180]);
    xlabel('Target location (°)'); ylabel('Final Saccade Error (MAE; °)');
    set(gca,'TickDir','out');
end


%% plot probability for R2choose against R1 precision as function of location bin  plus correlation 

my_Bin = 0:30:360; % create bin limit for 0 to 360, make 1 bin for 30 degree interval therefore result in 12 bins in total
%[BINS, EDGES] = discretize(X,N) where N is a scalar integer, divides the range of X into N uniform bins, and also returns the bin edges.
location_Bin = discretize(all_data.s_all.TarT360,my_Bin);%put target location data into bin from 1-12, e.g., 1 is 0-30 degree and so on.

bin_colors = hsv(length(my_Bin));%give each bin a distinct color

all_corr = nan(length(u_subj),1);
% figure();
% scatter(all_R1muErr(:,1),all_LOC_nChoose(:,1)); 
% H = lsline;
% [R,P] = corrcoef(all_R1muErr(:,1),all_LOC_nChoose(:,1));
% annotation('textbox',[.8 .8 .1 .1],'String',sprintf('R = %.02f, p = %.02f',R(2),P(2)))

figure;
for ss = 1:length(u_subj) % go through all the subjects

    this_n_choose = sum(all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(3) & all_data.use_trial==1);% compute all number of trials in choose condition

    all_R1muErr = nan(length(my_Bin)-1,1);
    all_LOC_nChoose = nan(length(my_Bin)-1,1);

    for bb = 1:(length(my_Bin)-1)



        subplot(4,5,ss); hold on;

        thisidx  = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(1) & all_data.use_trial==1 & location_Bin==bb;
        thisidx2 = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(3) & all_data.use_trial==1 & location_Bin==bb;

        this_x = mean(all_data.s_all.f_sacc_err(thisidx));% compute the mean of memory error for cue condition
        this_y = sum(thisidx2)/this_n_choose; % compute the proportion of choose trials in different location bin over the overall number of choose trial

        %scatter(this_x,this_y,'o','Color',bin_colors(bb,:),'MarkerFaceColor',bin_colors(bb,:));%,'MarkerFaceColor',cond_colors(cc,:));
        plot(this_x,this_y,'o','Color',bin_colors(bb,:),'MarkerFaceColor',bin_colors(bb,:),'MarkerSize',4);

        all_R1muErr(bb) = this_x;
        all_LOC_nChoose (bb) = this_y;

    end

    all_corr(ss) = corr(all_R1muErr,all_LOC_nChoose);

    text(4.5,0.27,u_subj{ss},'FontAngle','italic','HorizontalAlignment','right');


    xlim([0 5]);ylim([0 0.3]);
    xticks(0:2:6);yticks(0:0.1:0.3);

    % set(gca,'TickDir','out','XTick',[-180 0 180],'YTick',[-180 0 180]);
    if ss == 16; xlabel('R1 Error (MAE; °)'); ylabel('Proportion location chosen');end
    set(gca,'TickDir','out','LineWidth',0.75);

end

match_xlim(get(gcf,'Children')); match_ylim(get(gcf,'Children'));


all_corr_z = atanh(all_corr);

% compute mean and SEM using R-to-Z transformed values, plot these
% converted back to R (see below)
thismz = mean(all_corr_z);
thisez = std(all_corr_z)/sqrt(length(subj));


% plot all correlation values & their mean
figure; hold on;

plot([1.1;1.1],tanh(thismz+[-1;1]*thisez),'k-','LineWidth',1.5);
plot(ones(size(all_corr))-0.1,all_corr,'ko','MarkerFaceColor','k');
plot(1.1,tanh(thismz),'ko','MarkerSize',10,'LineWidth',1.5,'MarkerFaceColor','w');

xlim([0.5 1.5]); ylim([-1 1]);xticks(1); yticks(-1:0.5:1);
set(gca,'TickDir','out','LineWidth',0.75);
ylabel('Correlation');


% compute and spit out stats for this test
[~,thisp,~,thisstats] = ttest(all_corr_z);

fprintf('\nLocation heuristic analysis: t-test of correlations against 0\n');
fprintf('T-test: corr against zero, T(%i) = %.05f, p = %.05f, dz = %.05f\n\n',thisstats.df,thisstats.tstat,thisp,thisstats.tstat/sqrt(length(subj)));

clear thisp thisstats;


% plot a quick figure of location bins
figure; hold on;
for bb = 1:(length(my_Bin)-1)
    tmpth = linspace(my_Bin(bb),my_Bin(bb+1),101);

    plot(cosd(tmpth),sind(tmpth),'-','LineWidth',3,'Color',bin_colors(bb,:));
end
plot(0,0,'ko','MarkerFaceColor','k');
xlim([-1.1 1.1]); ylim([-1.1 1.1]); axis square; axis off; 


% figure();
% scatter(tmpErr,tmpP); 
% H = lsline;
% [R,P] = corrcoef(tmpErr,tmpP);
% annotation('textbox',[.8 .8 .1 .1],'String',sprintf('R = %.02f, p = %.02f',R(2),P(2)))
%% plot R2 Cue error against location choice see if correlation needed

all_R2muErr= nan(length(my_Bin)-1, length(u_subj));
figure;
for ss = 1:length(u_subj) % go through all the subjects

    this_n_choose = sum(all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(3) & all_data.use_trial==1);% compute all number of trials in choose condition


    for bb = 1:(length(my_Bin)-1)


        %subplot(length(my_Bin),length(u_subj),ss+length(u_subj)*(bb-1)); hold on;
        subplot(4,5,ss); hold on;

        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(2) & all_data.use_trial==1 & location_Bin==bb;
        thisidx2 = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(3) & all_data.use_trial==1 & location_Bin==bb;

        this_x = mean(all_data.s_all.f_sacc_err(thisidx));% compute the mean of memory error for cue condition
        this_y = sum(thisidx2)/this_n_choose;% compute the proportion of choose trials in different location bin over the overall number of choose trial

        plot(this_x,this_y,'o','Color',bin_colors(bb,:),'MarkerFaceColor',bin_colors(bb,:));
        %scatter(this_x,this_y,'o','Color',bin_colors(bb,:),'MarkerFaceColor',bin_colors(bb,:));%,'MarkerFaceColor',cond_colors(cc,:));
        all_LOC_nChoose(bb,ss) = this_y;
        all_R2muErr(bb,ss) = this_x;

    end

    text(4.5,0.35,u_subj{ss},'FontAngle','italic','HorizontalAlignment','right');


    xlim([0 5]);ylim([0 0.4]);
%    xticks(0:2:6);yticks(0:0.2:0.4);

    % set(gca,'TickDir','out','XTick',[-180 0 180],'YTick',[-180 0 180]);
    if ss == 16; xlabel('R1 Error (MAE; °)'); ylabel('Proportion of Choose');end
    set(gca,'TickDir','out','LineWidth',0.75);


    set(gca,'TickDir','out','XTick',[0:2:5],'YTick',[0:0.2:0.5]);
    if ss == 16; xlabel('R2 Cue Error (MAE; °)'); ylabel('Proportion location chosen');end
end

match_xlim(get(gcf,'Children')); match_ylim(get(gcf,'Children'));


%% plot error against rt see if correlation needed

all_corr_RT_err = nan(length(cu),length(u_subj));

figure;
for ss = 1:length(u_subj) % go through all the subjects

    %subplot(4,5,ss); hold on

    for cc = 1:length(cu)


        subplot(length(cu),length(u_subj),ss+length(u_subj)*(cc-1)); hold on;

        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(cc) & all_data.use_trial==1;
        thisidx2 = find(thisidx);

        %scatter(all_data.s_all.i_sacc_rt(thisidx2),all_data.s_all.f_sacc_err(thisidx2),'o','ColorVariable',);
        %plot(all_data.s_all.i_sacc_rt(thisidx2),all_data.s_all.f_sacc_err(thisidx2),'o','Color',cond_colors(cc,:),'MarkerFaceColor',cond_colors(cc,:));
        scatter(all_data.s_all.i_sacc_rt(thisidx2),all_data.s_all.f_sacc_err(thisidx2),10,cond_colors(cc,:),'filled','MarkerFaceAlpha',0.3);%,'MarkerFaceColor',cond_colors(cc,:));
        [R,P] = corrcoef(all_data.s_all.i_sacc_rt(thisidx2),all_data.s_all.f_sacc_err(thisidx2));

        all_corr_RT_err(cc,ss) = R(1,2);
     
       
    end


    xlim([0 1.5]);ylim([0 5]);
    set(gca,'TickDir','out','XTick',[0 .5 1.5],'YTick',[0 2.5 5]);
    xlabel('Response Time'); ylabel('Final Saccade Error');
end

match_xlim(get(gcf,'Children')); match_ylim(get(gcf,'Children'));


% Accuracy & RT correlation bin by conditions
figure;
for ss = 1:length(u_subj) % go through all the subjects

    %subplot(4,5,ss); hold on
    subplot(1,1,1);hold on
    for cc = 1:length(cu)


        thisidx = all_data.subj_all==ss & all_data.s_all.trialinfo(:,1)==cu(1) & all_data.use_trial==1;
        thisidx2 = find(thisidx);

        %scatter(all_data.s_all.i_sacc_rt(thisidx2),all_data.s_all.f_sacc_err(thisidx2),'o','ColorVariable','MarkerFaceAlpha',0.3);
        %plot(all_data.s_all.i_sacc_rt(thisidx2),all_data.s_all.f_sacc_err(thisidx2),10,'Color',cond_colors(cc,:),'MarkerFaceAlpha',0.3);
        scatter(all_data.s_all.i_sacc_rt(thisidx2),all_data.s_all.f_sacc_err(thisidx2),10,cond_colors(1,:),'Filled','MarkerFaceColor',cond_colors(1,:));%,'MarkerFaceColor',cond_colors(cc,:));
        %H = lsline;
        %[R,P] = corrcoef(all_data.s_all.i_sacc_rt(thisidx2),all_data.s_all.f_sacc_err(thisidx2));
        %annotation('textbox',[.8 .8 .1 .1],'String',sprintf('R = %.02f, p = %.02f',R(2),P(2)))


        xlim([0 1.5]);ylim([0 5]);
        set(gca,'TickDir','out','XTick',[0 .5 1.5],'YTick',[0 2.5 5]);
        xlabel('Response Time (s)'); ylabel('Final Saccade Error');

    end

end

[R,P] = corrcoef(all_data.s_all.i_sacc_rt(thisidx2),all_data.s_all.f_sacc_err(thisidx2));
annotation('textbox',[.8 .8 .1 .1],'String',sprintf('R = %.02f, p = %.02f',R(2),P(2)))

match_xlim(get(gcf,'Children')); match_ylim(get(gcf,'Children'));

% Accuracy & RT individual correlation bin by conditions 

SEM_1  = std(all_corr_RT_err(1,:))/sqrt(length(all_corr_RT_err(1,:)));
SEM_2  = std(all_corr_RT_err(2,:))/sqrt(length(all_corr_RT_err(2,:)));
SEM_3  = std(all_corr_RT_err(3,:))/sqrt(length(all_corr_RT_err(3,:)));

figure;
for ss = 1:length(all_corr_RT_err)

    subplot(1,1,1);hold on
    scatter(1,all_corr_RT_err(1,:), 'o','MarkerFaceColor',cond_colors(1,:),'MarkerEdgeColor',cond_colors(1,:),'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
    scatter(2,all_corr_RT_err(2,:), 'o','MarkerFaceColor',cond_colors(2,:),'MarkerEdgeColor',cond_colors(2,:),'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
    scatter(3,all_corr_RT_err(3,:), 'o','MarkerFaceColor',cond_colors(3,:),'MarkerEdgeColor',cond_colors(3,:),'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
 
    errorbar(1,mean(all_corr_RT_err(1,:)),SEM_1,'--ko','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black','CapSize',10);
    errorbar(2,mean(all_corr_RT_err(2,:)),SEM_2,'--ko','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black','CapSize',10);
    errorbar(3,mean(all_corr_RT_err(3,:)),SEM_3,'--ko','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black','CapSize',10);

   
    set(gca,'XTick',1:length(cu),'TickDir','out','LineWidth',1.5,'XTickLabel',cond_str,'FontSize',14,'XTickLabelRotation',-45);
    xlim([0.5 0.5+length(cu)]);ylim([-.5 .5]);ylabel('Correlation coefficient');
    plot([0 4],[0 0],'k--','LineWidth',1.5);
end
    
 match_xlim(get(gcf,'Children')); match_ylim(get(gcf,'Children'));

fprintf('\nRT vs error for each condition\n');

for cc = 1:size(all_corr_RT_err,1)
    all_corr_RT_err_z = atanh(all_corr_RT_err(cc,:));

    [~,thisp,~,thisstats] = ttest(all_corr_RT_err_z);

    fprintf('T-test: %s against 0, T(%i) = %.03f, p = %0.05f, d = %0.05f\n',cond_str{cc}, thisstats.df,thisstats.tstat,thisp,mean(all_corr_RT_err_z)/std(all_corr_RT_err_z));



end 
%% stats - shuffle condition labels within each subj before computing distributions, F-scores
% - use only included trials? yes

fprintf('\nStats starting!!\n');



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
            
%             if ii == 1
%                 fprintf('\n\nRunning stats for %s - %s\n',to_plot_2d{pp},dim_str{dim_to_plot(dd)});
%             end

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

            % do a full 'real' stats w/ anovan for the first iteration (to
            % match euclidean error at top)
            if ii == 1
                % stats for this parameter (ANOVA)
                [thisp,thisanovatab] = anovan( thisX(:,1),thisX(:,[2 3])  , 'random',2,'varnames',{'condition','subj'},'model','interaction' ,'display','off');

                % compute partial eta squared (SS_cond/(SS_cond+SS_err))
                this_peta2 = thisanovatab{2,2}/(thisanovatab{2,2}+thisanovatab{4,2}); % SS_cond/(SS_cond + SS_cond*subj)

                fprintf('\nParametric stats - %s, %s\n',to_plot_2d{pp},dim_str{dim_to_plot(dd)});
                fprintf('1-way RM AOV: F(%i,%i) = %.03f, p = %0.05f, partial eta^2 = %0.05f\n',thisanovatab{2,3},thisanovatab{4,3},thisanovatab{2,6},thisp(1),this_peta2);
            end

            
            for cp_idx = 1:size(cond_pairs,1)
%                 if ii == 1
%                     allT{cp_idx} = nan(length(to_plot_2d),length(dim_to_plot),niter+1);
%                 end
                [~,tmp_p,~,tmp_stats] = ttest(thisX(thisX(:,2)==cond_pairs(cp_idx,1),1),thisX(thisX(:,2)==cond_pairs(cp_idx,2),1));
                if ii == 1
                    fprintf('%s vs %s: T(%i) = %0.03f, p = %0.03f, dz = %0.05f\n',cond_str{cond_pairs(cp_idx,1)},cond_str{cond_pairs(cp_idx,2)},tmp_stats.df,tmp_stats.tstat,tmp_p,tmp_stats.tstat/sqrt(length(subj)));
                end
                allT{cp_idx}(pp,dd,ii) = tmp_stats.tstat;
                clear tmp_stats tmp_p;
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

fprintf('\n\n');

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

            % do a full 'real' stats w/ anovan for the first iteration (to
            % match euclidean error at top)
            if ii == 1
                % stats for this parameter (ANOVA)
                [thisp,thisanovatab] = anovan( thisX(:,1),thisX(:,[2 3])  , 'random',2,'varnames',{'condition','subj'},'model','interaction' ,'display','off');

                % compute partial eta squared (SS_cond/(SS_cond+SS_err))
                this_peta2 = thisanovatab{2,2}/(thisanovatab{2,2}+thisanovatab{4,2}); % SS_cond/(SS_cond + SS_cond*subj)

                fprintf('\nParametric stats - %s, avg\n',to_plot_2d{pp});
                fprintf('1-way RM AOV: F(%i,%i) = %.03f, p = %0.05f, partial eta^2 = %0.05f\n',thisanovatab{2,3},thisanovatab{4,3},thisanovatab{2,6},thisp(1),this_peta2);
            end
            
            for cp_idx = 1:size(cond_pairs,1)
%                 if ii == 1
%                     allT{cp_idx} = nan(length(to_plot_2d),length(dim_to_plot),niter+1);
%                 end
                [~,tmp_p,~,tmp_stats] = ttest(thisX(thisX(:,2)==cond_pairs(cp_idx,1),1),thisX(thisX(:,2)==cond_pairs(cp_idx,2),1));
                if ii == 1
                    fprintf('%s vs %s: T(%i) = %0.03f, p = %0.03f, dz = %0.05f\n',cond_str{cond_pairs(cp_idx,1)},cond_str{cond_pairs(cp_idx,2)},tmp_stats.df,tmp_stats.tstat,tmp_p,tmp_stats.tstat/sqrt(length(subj)));
                end

                allT{cp_idx}(pp,ii) = tmp_stats.tstat;
                clear tmp_stats tmp_p;
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