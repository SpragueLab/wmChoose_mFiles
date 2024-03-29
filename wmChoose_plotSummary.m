% wmChoose_plotSummary
%
% plots summary of data from each subj of wmChoose experiment, does not
% sort by condition
%
% for each subj, several figures:
% - plot of eye position timecourse with targets overlaid (2d) [DONE]
% - plot of target presentation; delay timecourse (indicate whether trial
%   excluded) [OR: image of this for X, Y, absolute(or signed) diff from
%   0?) scaled to value used for cutoff?
% - overall summary by run of type of trial exclusion for each trial
% - histogram of error for initial, final saccade

% EXCLUSION LABELS:
% - 11: drift correction too big [d]
% - 12: calibration out of range (can't visualize either of these here..)[c]
% - 13: fixation outside of range during delay epoch [fix red]
% - 20: no primary saccade found [fix x]
% - 21: duration too long/amplitude too small [red traj]
% - 22: error too large (primary) [red traj]

root = '/Volumes/data/wmChoose';

%subj = {'aa1','aa2','ab1','ab2','ac1','ac2','ae','af','ag'};
subj = {'ai'};

QC_dir = 'preproc_QC';

% which trials do we want to exclude in the final sumamry plot?
% (for report, etc, will still use all criteria, at least for each one
% individually, but we may not wish to actually drop trials based on
% these)
WHICH_EXCL = [11 13 20 21 22]; % NOTE: will always show all exclusion criteria on each subplot, but

% for final exclusion report at the end
all_excl = [11 12 13 20 21 22];
excl_labels = {'drift','calibration','delay fixation','no i_sacc','bad i_sacc','i_sacc err'};

SAMPLING_RATE = 1000;
SAMPLING_PER  = 1/SAMPLING_RATE;


% for plotting: how many rows/cols for 2d trial plots (hopefully a bunch!)
NROWS_2D = 10;
NROWS_1D = 30;
NCOLS_1D = 10;

NROWS_FULLTRIAL = 8;
NCOLS_FULLTRIAL = 10; % for plotting full trial timecourse

tmp_lines = lines(7);

TARG_COLORS = lines(5);
TARG_COLORS = TARG_COLORS(4:5,:);
FIG_POS_2D = [34 501 2508 816];
FIG_POS_1D = [1000 120 1402 1218];
SACC_COLORS = [0 0 0; tmp_lines(3,:)];

EXCL_COLOR = [0.8 0 0];

% XDAT tags that are useful
XDAT_RESP = 3;
XDAT_PRERESP = [1 2];
XDAT_ITI = 5;

RAW_COLORS = lines(2); % for 1d traces (X,Y)

for ss = 1:length(subj)
    
    fn = sprintf('%s/data/%s_wmChoose_behav.mat',root,subj{ss});
    fprintf('loading %s\n',fn);
    this_data = load(fn);
    clear fn;
    
    ru = unique(this_data.r_all);
    nfigs = ceil(length(ru)/NROWS_2D);
    
    % 2d plot of each trial (for now, just trace & coords)
    % - also include hints as to why the trial was excluded (fix point, or
    %   trace will be red; also text markers - see comments below & key
    %   above)
    figs_2d(1) = figure; set(gcf,'Position',FIG_POS_2D); fig_cnt = 1; %row_cnt = 1; col_cnt = 1; fig_cnt =1 ;
    set(gcf,'NumberTitle','off','Name',sprintf('%s (%i of %i)',subj{ss},fig_cnt,nfigs));
    for ii = 1:size(this_data.c_all,1)
        
        thisrow = mod(find(ru==this_data.r_all(ii)),NROWS_2D);
        if thisrow==0
            thisrow = NROWS_2D; % for subplot selection purposes...
        end
        
        if ii == 1
            NCOLS_2D = max(this_data.t_all);
        end
        thiscol = this_data.t_all(ii);
        
        %         if thiscol==1
        %             thisrow
        %             this_data.r_all(ii)
        %         end
        
        
        subplot(NROWS_2D,NCOLS_2D,(thisrow-1)*NCOLS_2D+thiscol); hold on;
        
        % if trial is excluded (per our chosen criteria!), mark with red ring; otherwise, black
        if any(ismember(this_data.s_all.excl_trial{ii},WHICH_EXCL))
            plot(15*cos(linspace(2*pi/180,2*pi,180)), 15*sin(linspace(2*pi/180,2*pi,180)),'-','Color',EXCL_COLOR);
        else
            plot(15*cos(linspace(2*pi/180,2*pi,180)), 15*sin(linspace(2*pi/180,2*pi,180)),'-','Color',[0 0 0]);
        end
        
        % draw both targets
        for targ_idx = 1:2
            plot(this_data.coords_all{targ_idx}(ii,1),this_data.coords_all{targ_idx}(ii,2),'o','MarkerSize',5,'Color',TARG_COLORS(targ_idx,:),'MarkerFaceColor',TARG_COLORS(targ_idx,:));
        end
        
        % if there's a saccade to draw, draw it
        if ~isempty(this_data.s_all.i_sacc_trace{ii})
            
            % and if the trial were excluded for a saccade-related reason,
            % mark it w/ red i_sacc trace
            if any(ismember(this_data.s_all.excl_trial{ii},[21 22]))
                plot(this_data.s_all.i_sacc_trace{ii}(:,1),this_data.s_all.i_sacc_trace{ii}(:,2),'-','LineWidth',1.5,'Color',EXCL_COLOR);
                % otherwise, black
            else
                plot(this_data.s_all.i_sacc_trace{ii}(:,1),this_data.s_all.i_sacc_trace{ii}(:,2),'-','LineWidth',1.5,'Color',SACC_COLORS(1,:));
            end
        end
        
        % make sure there's a unique final saccade (initial saccade also
        % included as final when n_sacc == 1)
        if ~isempty(this_data.s_all.f_sacc_trace{ii}) && this_data.s_all.n_sacc(ii) > 1 % make sure we don't overlay initial/final
            plot(this_data.s_all.f_sacc_trace{ii}(:,1),this_data.s_all.f_sacc_trace{ii}(:,2),'-','LineWidth',1.5,'Color',SACC_COLORS(2,:));
        end
        
        % by default, we'll mark fix w/ a +, but use an x if excluded for
        % absence of primary saccade detected
        fixmarker = '+';
        
        % mark fix red if no primary saccade detected or fixation break
        % during delay
        if any(ismember(this_data.s_all.excl_trial{ii},[13 20]))
            fixcolor = EXCL_COLOR;
            if any(ismember(this_data.s_all.excl_trial{ii},20))
                fixmarker = 'x';
            end
        else
            fixcolor = [0 0 0];
        end
        plot(0,0,fixmarker,'Color',fixcolor);
        clear fixmarker fixcolor;
        
        % add a marker for preproc-related reasons for dropping a trial
        % (failed calibration/drift correction according to thresholds)
        exclstr = '';
        if any(ismember(this_data.s_all.excl_trial{ii},11))
            exclstr(end+1) = 'd';
        end
        if any(ismember(this_data.s_all.excl_trial{ii},12))
            exclstr(end+1) = 'c';
        end
        if ~isempty(exclstr)
            text(-14.7,14.7,exclstr,'Color',EXCL_COLOR,'VerticalAlignment','middle');
        end
        clear exclstr
        
        axis square off;
        xlim([-15 15]);ylim([-15 15]);
        
        
        if thisrow==NROWS_2D && thiscol==NCOLS_2D && ii~=size(this_data.c_all,1)
            figs_2d(end+1) = figure; set(gcf,'Position',FIG_POS_2D); fig_cnt = fig_cnt+1; % we'll need a new figure for next iteration
            set(gcf,'NumberTitle','off','Name',sprintf('%s (%i of %i)',subj{ss},fig_cnt,nfigs));
        end
        
        clear thiscol thisrow
        
    end
    clear NCOLS_2D;
    
    for ff = 1:length(figs_2d)
        saveas(figs_2d(ff),sprintf('%s/preproc_QC/%s_eachTrial_2d_%i_of_%i.png',root,subj{ss},ff,nfigs));
    end
    clear figs_2d;
    

    
    % let's also try it the imagesc way...
    figs_1d(1) = figure;set(gcf,'Position',FIG_POS_1D);
    for ii = 1:size(this_data.c_all,1)
        
        thisrow = 1;
        thiscol = this_data.r_all(ii);
        
        thiscol = mod(thiscol,NCOLS_1D);
        if thiscol==0
            thiscol = NCOLS_1D;
        end
        
        if thiscol==1 && this_data.r_all(ii)>NCOLS_1D && this_data.t_all(ii)==1
            % if this is the first trial of a new run on a new figure
            figs_1d(end+1)=figure;set(gcf,'Position',FIG_POS_1D);
        end
        
        
        subplot(1,NCOLS_1D,(thisrow-1)*NCOLS_1D+thiscol); hold on;
        
        preresp_idx = ismember(this_data.s_all.XDAT{ii},XDAT_PRERESP);
        
        
        imagesc(0.001*(1:sum(preresp_idx)),this_data.t_all(ii),sqrt( this_data.s_all.X{ii}(preresp_idx).^2 + this_data.s_all.Y{ii}(preresp_idx).^2 ).');
        
        if any(ismember(this_data.s_all.excl_trial{ii},WHICH_EXCL)) % if we do exclude
            % draw an x to the right of the trial
            plot(4.25,this_data.t_all(ii),'x','LineWidth',1.5,'Color',EXCL_COLOR)
            if any(ismember(this_data.s_all.excl_trial{ii},13))
                % if delay out-of-bounds, add a circle
                plot(4.25,this_data.t_all(ii),'o','LineWidth',1.5,'Color',EXCL_COLOR)
            end
        end
        
        set(gca,'CLim',[0 1]*this_data.excl_criteria.delay_fix_thresh);
        
        axis ij off tight;
    end
    for ff = 1:length(figs_1d)
        saveas(figs_1d(ff),sprintf('%s/preproc_QC/%s_delay_%i_of_%i.png',root,subj{ss},ff,length(figs_1d)));
    end
    clear figs_1d;
    
    % now let's overlay all aligned traces (subplots for included/excluded)
    subplot_group{1} = [1 3 5];
    subplot_group{2} = [2 4 6];
    figure; subplot(4,2,subplot_group{1}); hold on;
    plot(0,0,'ks','MarkerSize',8,'MarkerFaceColor','k');
    subplot(4,2,subplot_group{2}); hold on;
    plot(0,0,'ks','MarkerSize',8,'MarkerFaceColor','k');
    
    for ii = 1:size(this_data.c_all,1)
        if any(ismember(this_data.s_all.excl_trial{ii},WHICH_EXCL))
            thiscolor = EXCL_COLOR; subplot(4,2,subplot_group{2});
        else
            thiscolor = [0.2 0.2 0.2]; subplot(4,2,subplot_group{1});
        end
        
        
        if ~isempty(this_data.s_all.i_sacc_trace{ii}) && ~isempty(this_data.s_all.f_sacc_trace{ii})
            % initial saccade
            [tmpth,tmpr] = cart2pol(this_data.s_all.i_sacc_trace{ii}(:,1), this_data.s_all.i_sacc_trace{ii}(:,2));
            
            % for cued response trials, always rotate to cued stimulus
            if this_data.c_all(ii,1)~=3
                thisoffset = this_data.c_all(ii,2);
            else
                % if uncued, rotate to the 'selected' target
                thisoffset = this_data.c_all(ii,1+this_data.s_all.sel_targ(ii));
            end
            
            tmpth = tmpth-deg2rad(thisoffset);
            [alignx,aligny] = pol2cart(tmpth,tmpr);
            plot(alignx,aligny,'-','Color',thiscolor);
            clear alignx aligny tmpth tmpr;
            
            if this_data.s_all.n_sacc(ii)>1
                [tmpth,tmpr] = cart2pol(this_data.s_all.f_sacc_trace{ii}(:,1), this_data.s_all.f_sacc_trace{ii}(:,2));
                tmpth = tmpth-deg2rad(thisoffset);
                [alignx,aligny] = pol2cart(tmpth,tmpr);
                plot(alignx,aligny,'-','Color',thiscolor);
                clear alignx aligny tmpth tmpr;
            end
            
        end
        
        if any(ismember(this_data.s_all.excl_trial{ii},WHICH_EXCL))
            subplot(4,2,8); hold on;
        else
            subplot(4,2,7); hold on;
        end
        % also draw the full timeseries of r
        [~,tmpr] = cart2pol(this_data.s_all.X{ii},this_data.s_all.Y{ii});
        thist = (0.001)*(1:length(tmpr)); % 1000 Hz...
        startt = thist(find(this_data.s_all.XDAT{ii}==XDAT_RESP,1,'first'))-0.5;
        drawidx = this_data.s_all.XDAT{ii}~=XDAT_ITI & thist.'>=startt;
        plot(thist(drawidx)-startt-0.5,tmpr(drawidx),'-','Color',thiscolor);
        %TODO: also draw the saccades on top of this?
        
        
    end
    for ii = 1:2
        subplot(4,2,subplot_group{ii}); % finish off each subplot
        plot(15*cosd(linspace(360/1000,360,1000)),15*sind(linspace(360/1000,360,1000)),'k-')
        plot(12,0,'o','Color',[0 0 0.7],'MarkerSize',10,'MarkerFaceColor',[0 0 0.7]);
        xlim([-15 15]); ylim([-15 15]); axis square off;
        subplot(4,2,6+ii);
        xlim([-0.5 2.5]);
        ylim([-3 15]);
        plot([-0.5 2.5],[0 0],'k-');
        plot([-0.5 2.5],[12 12],'k--');
        %plot(-0.4*[1 1],[0 1]+5,'k-','LineWidth',2)
        set(gca,'XTick',0:2,'YTick',[0 12],'TickDir','out');
    end
    saveas(gcf,sprintf('%s/preproc_QC/%s_allTraces.png',root,subj{ss}));
    
    
    
    %% ~~~~~ plot full trial timecourse, highlighting primary/final saccade
    
    % maybe do 10 trials across, 8 tall?
    % - first, plot full X, Y timecourse
    % - then,  use latency to plot X_trace, Y_trace of init, final saccade
    
    nfigs = ceil(size(this_data.c_all,1)/(NCOLS_FULLTRIAL*NROWS_FULLTRIAL));
    figs_ft(1) = figure; set(gcf,'Position',FIG_POS_2D); fig_cnt = 1; %row_cnt = 1; col_cnt = 1; fig_cnt =1 ;
    set(gcf,'NumberTitle','off','Name',sprintf('Full trial - %s (%i of %i)',subj{ss},fig_cnt,nfigs));
    for ii = 1:size(this_data.c_all,1)
        
        thisrow = mod(ceil(ii/NCOLS_FULLTRIAL),NROWS_FULLTRIAL);
        if thisrow==0
            thisrow = NROWS_FULLTRIAL; % for subplot selection purposes...
        end
        
        thiscol = mod(ii,NCOLS_FULLTRIAL);
        if thiscol == 0
            thiscol = NCOLS_FULLTRIAL;
        end
        %         if thiscol==1
        %             thisrow
        %             this_data.r_map(ii)
        %         end
        
        
        subplot(NROWS_FULLTRIAL,NCOLS_FULLTRIAL,(thisrow-1)*NCOLS_FULLTRIAL+thiscol); hold on;
        
        myt = (1:length(this_data.s_all.X{ii})) * SAMPLING_PER;
        
        
        % when did response epoch start? (XDAT=4)
        this_resp_start = myt(find(this_data.s_all.XDAT{ii}==XDAT_RESP,1,'first'));
        this_resp_end   = myt(find(this_data.s_all.XDAT{ii}==XDAT_RESP,1,'last'));
        
        
        % start by just plotting raw X, Y for this trial
        
        % TARG positions
        plot([myt(1) myt(end)],[1 1]*this_data.coords_all{1}(ii,1),'--','Color',TARG_COLORS(1,:));
        plot([myt(1) myt(end)],[1 1]*this_data.coords_all{1}(ii,2),'--','Color',TARG_COLORS(1,:));
        plot([myt(1) myt(end)],[1 1]*this_data.coords_all{2}(ii,1),'--','Color',TARG_COLORS(2,:));
        plot([myt(1) myt(end)],[1 1]*this_data.coords_all{2}(ii,2),'--','Color',TARG_COLORS(2,:));

        plot([myt(1) myt(end)],[0 0],'k-');
        
        % plot start of response epoch
        plot(this_resp_start*[1 1],[-15 15],'-','Color',[0.5 0.5 0.5]);
        plot(this_resp_end*[1 1],  [-15 15],'-','Color',[0.5 0.5 0.5]);
        
        
        % X, Y eye positions
        plot(myt,this_data.s_all.X{ii},'-','LineWidth',1,'Color',RAW_COLORS(1,:));
        plot(myt,this_data.s_all.Y{ii},'-','LineWidth',1,'Color',RAW_COLORS(2,:));
        
        
        % overlay the primary saccade
        i_sacc_t = (1:size(this_data.s_all.i_sacc_trace{ii},1)) * SAMPLING_PER + this_resp_start + this_data.s_all.i_sacc_rt(ii);
        if ~isempty(this_data.s_all.i_sacc_trace{ii})
            plot(i_sacc_t,this_data.s_all.i_sacc_trace{ii}(:,1),'-','LineWidth',2,'Color',SACC_COLORS(1,:));
            plot(i_sacc_t,this_data.s_all.i_sacc_trace{ii}(:,2),'-','LineWidth',2,'Color',SACC_COLORS(1,:));
        end
        
        if this_data.s_all.n_sacc(ii) >= 2
            f_sacc_t = (1:size(this_data.s_all.f_sacc_trace{ii},1)) * SAMPLING_PER + this_resp_start + this_data.s_all.f_sacc_rt(ii);
            if ~isempty(this_data.s_all.f_sacc_trace{ii})
                plot(f_sacc_t,this_data.s_all.f_sacc_trace{ii}(:,1),'-','LineWidth',1.5,'Color',SACC_COLORS(2,:));
                plot(f_sacc_t,this_data.s_all.f_sacc_trace{ii}(:,2),'-','LineWidth',1.5,'Color',SACC_COLORS(2,:));
            end
            
        end
        % run, trial, sess # (red if excluded)
        % if actually exlcude (if ismember this trial, which_excl), red
        % if satisfies an exclusion criterion (even if we don't exclude),
        % italics
        this_angle = 'normal';
        this_color = [0 0 0];
        if ~isempty(this_data.s_all.excl_trial{ii})
            this_angle = 'italic';
            if any(ismember(this_data.s_all.excl_trial{ii},WHICH_EXCL))
                this_color = EXCL_COLOR;
            end
        end
        text(0,15,sprintf('r%02.f, %02.f',this_data.r_all(ii),this_data.t_all(ii)),'FontSize',10,'FontAngle',this_angle,'Color',this_color);
        
        set(gca,'XLim',[0 15],'YLim',[-15 15]);
        axis off;
        
        
        if thisrow==NROWS_FULLTRIAL && thiscol==NCOLS_FULLTRIAL && ii~=size(this_data.c_all,1)
            figs_ft(end+1) = figure; set(gcf,'Position',FIG_POS_2D); fig_cnt = fig_cnt+1; % we'll need a new figure for next iteration
            set(gcf,'NumberTitle','off','Name',sprintf('Full trial - %s (%i of %i)',subj{ss},fig_cnt,nfigs));
        end
        
        clear thiscol thisrow
        
    end
    for ff = 1:length(figs_ft)
        saveas(figs_ft(ff),sprintf('%s/%s/%s_fullTrial_%i_of_%i.png',root,QC_dir,subj{ss},ff,length(figs_ft)));
    end
    clear figs_ft;
    
    
    
    
    
    
    
    %%
    % exclusion report (bar graph; imgs)
    % first, make a subplot for each run and show n_trials x n_excl dots
    % indicating why a given trial was dropped (color it red if we use that
    % criterion, black if it faield, but we're ignoring it)
    figure;
    ru = unique(this_data.r_all);
    NCOLS_EXCL = 10;
    NROWS_EXCL = ceil(length(ru)/NCOLS_EXCL);
    
    for rr = 1:length(ru)
        
        subplot(NROWS_EXCL,NCOLS_EXCL,rr); hold on; axis ij;
        
        % which trials in this run?
        runidx = find(this_data.r_all==ru(rr));
        
        
        for tt = 1:length(runidx)
            % if there's something to plot...
            if ~isempty(this_data.s_all.excl_trial{runidx(tt)})
                
                thisexcl = this_data.s_all.excl_trial{runidx(tt)};
                thisy = this_data.t_all(runidx(tt));
                
                for ee = 1:length(thisexcl)
                    
                    % x value is position in all_excl
                    thisx = find(all_excl==thisexcl(ee));
                    
                    % color is whether value is in WHICH_EXCL
                    if ismember(thisexcl(ee),WHICH_EXCL)
                        thiscolor = EXCL_COLOR;
                    else
                        thiscolor = [0 0 0];
                    end
                    
                    plot(thisx,thisy,'o','MarkerSize',5,'MarkerFaceColor',thiscolor,'Color',thiscolor);
                    
                end
                clear thisx thisy thiscolor;
            end
            
        end
        
        set(gca,'YTick',[1 5:5:length(runidx)]);
        
        title(sprintf('Run %i',ru(rr)));
        if rr==1
            ylabel('Trial');
            set(gca,'YTickLabel',[1 5:5:length(runidx)]);
        else
            set(gca,'YTickLabel',[]);
        end
        
        set(gca,'XTick',1:length(all_excl),'XTickLabel',excl_labels','XTickLabelRotation',-90,'TickDir','out','TickLabelInterpreter','none');
        xlim([0 length(all_excl)+1]);
        ylim([0 length(runidx)+1]);
        
        clear runidx;
        
    end
    set(gcf,'Position',[1000         881        1368         457]);
    saveas(gcf,sprintf('%s/preproc_QC/%s_excludedTrials.png',root,subj{ss}));
    
    % now just a quick bar graph across all runs of % trials falling under
    % each exclusion criterion (note, of course, that there will be
    % overlap)
    
    % TODO:
    figure;
    % two subplots: one of each criterion, color-coded as to whether it's
    % used
    subplot(1,2,1); hold on;
    for ee = 1:length(all_excl)
        if ismember(all_excl(ee),WHICH_EXCL)
            thiscolor = EXCL_COLOR;
        else
            thiscolor = [0 0 0];
        end
        
        this_excl = cellfun(@any,cellfun(@(a) ismember(a,all_excl(ee)),this_data.s_all.excl_trial,'UniformOutput',false));
        bar(ee,100*mean(this_excl),'FaceColor',thiscolor);
        clear this_excl;
    end
    set(gca,'XTick',1:length(all_excl),'XTickLabel',excl_labels,'XTickLabelRotation',-90);
    ylim([0 100]);
    
    % second, 2 bars, one for applying all criterion (black) and chosen
    % criterion (white)
    subplot(1,2,2); hold on;
    this_excl = ~cellfun(@isempty,this_data.s_all.excl_trial);
    bar(1,100*mean(this_excl),'EdgeColor',[0 0 0],'FaceColor',[0 0 0],'LineWidth',1.5);
    clear this_excl;
    
    this_excl = cellfun(@any,cellfun(@(a) ismember(a,WHICH_EXCL),this_data.s_all.excl_trial,'UniformOutput',false));
    bar(2,100*mean(this_excl),'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineWidth',1.5);
    set(gca,'XTick',[1 2],'XTickLabel',{'All','Specified'},'XTickLabelRotation',-45);
    xlabel('Exclusions');
    
    saveas(gcf,sprintf('%s/preproc_QC/%s_proportionExcluded.png',root,subj{ss}));
    
    
    
    clear this_data;
    
    if length(subj)>1
        close all;
    end
    
end