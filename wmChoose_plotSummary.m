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

root = '/Volumes/data/wmChoose';

subj = {'ae'};%,'aa2','ab1','ab2','ac1','ac2','ae','af','ag'}; 

% for plotting: how many rows/cols for 2d trial plots (hopefully a bunch!)
NROWS_2D = 10;
NROWS_1D = 30; 
NCOLS_1D = 10;

TARG_COLORS = lines(5);
TARG_COLORS = TARG_COLORS(4:5,:);
FIG_POS_2D = [34 501 2508 816];
FIG_POS_1D = [1000 120 1402 1218];
SACC_COLORS = [0 0 0; lines(1)];

RAW_COLORS = lines(2); % for 1d traces (X,Y)

for ss = 1:length(subj)
    
    fn = sprintf('%s/data/%s_wmChoose_behav.mat',root,subj{ss});
    fprintf('loading %s\n',fn);
    this_data = load(fn);
    clear fn;
    
    ru = unique(this_data.r_all);
    nfigs = ceil(length(ru)/NROWS_2D);
    
    % 2d plot of each trial (for now, just trace & coords)
    figure; set(gcf,'Position',FIG_POS_2D); fig_cnt = 1; %row_cnt = 1; col_cnt = 1; fig_cnt =1 ;
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
        
        subplot(NROWS_2D,NCOLS_2D,(thisrow-1)*NCOLS_2D+thiscol); hold on;
        if ~isempty(this_data.s_all.excl_trial{ii}) % if exclude
            plot(0,0,'r+');
        else % other wise
            plot(0,0,'k+');
        end
        plot(15*cos(linspace(2*pi/180,2*pi,180)), 15*sin(linspace(2*pi/180,2*pi,180)),'k-');
        for targ_idx = 1:2
            plot(this_data.coords_all{targ_idx}(ii,1),this_data.coords_all{targ_idx}(ii,2),'o','MarkerSize',5,'Color',TARG_COLORS(targ_idx,:),'MarkerFaceColor',TARG_COLORS(targ_idx,:));
        end
        if ~isempty(this_data.s_all.i_sacc_trace{ii})
            plot(this_data.s_all.i_sacc_trace{ii}(:,1),this_data.s_all.i_sacc_trace{ii}(:,2),'-','LineWidth',1.5,'Color',SACC_COLORS(1,:));
        end
        if ~isempty(this_data.s_all.f_sacc_trace{ii}) && this_data.s_all.n_sacc(ii) > 1 % make sure we don't overlay initial/final
            plot(this_data.s_all.f_sacc_trace{ii}(:,1),this_data.s_all.f_sacc_trace{ii}(:,2),'-','LineWidth',1.5,'Color',SACC_COLORS(2,:));
        end
        
        axis square off;
        xlim([-15 15]);ylim([-15 15]);
       

        if thisrow==NROWS_2D && thiscol==NCOLS_2D && ii~=size(this_data.c_all,1)
            figure; set(gcf,'Position',FIG_POS_2D); fig_cnt = fig_cnt+1; % we'll need a new figure for next iteration
            set(gcf,'NumberTitle','off','Name',sprintf('%s (%i of %i)',subj{ss},fig_cnt,nfigs));        
        end

        clear thiscol thisrow
        
    end
    clear NCOLS_2D;
    
    % 1d plot of each trial, pre-go-cue, colored accoridng to whether trial
    % excluded
    
    figure; set(gcf,'Position',FIG_POS_1D);fig_cnt = 1;
    for ii = 1:size(this_data.c_all,1)
        thisrow = this_data.t_all(ii);
        thiscol = this_data.r_all(ii);

        thiscol = mod(thiscol,NCOLS_1D);
        if thiscol==0
            thiscol = NCOLS_1D;
        end
        
        subplot(NROWS_1D,NCOLS_1D,(thisrow-1)*NCOLS_1D+thiscol); hold on;
        
        preresp_idx = ismember(this_data.s_all.XDAT{ii},[1 2]);
        
        plot(0.001*(1:sum(preresp_idx)).',zeros(sum(preresp_idx),1),'k--');
        
        if isempty(this_data.s_all.excl_trial{ii})
            % draw solid
            plot(0.001*(1:sum(preresp_idx)).',this_data.s_all.X{ii}(preresp_idx),'-','LineWidth',1.5,'Color',RAW_COLORS(1,:));
            plot(0.001*(1:sum(preresp_idx)).',this_data.s_all.Y{ii}(preresp_idx),'-','LineWidth',1.5,'Color',RAW_COLORS(2,:));
        else
            % draw dashed?
            plot(0.001*(1:sum(preresp_idx)).',this_data.s_all.X{ii}(preresp_idx),':','LineWidth',1.5,'Color',RAW_COLORS(1,:));
            plot(0.001*(1:sum(preresp_idx)).',this_data.s_all.Y{ii}(preresp_idx),':','LineWidth',1.5,'Color',RAW_COLORS(2,:));
        end
        
        ylim([-3 3]); xlim([0 4]);
        axis off;
        
        % if we're on the last row, col and not last trial, open a new fig
        if thisrow==NROWS_1D && thiscol==NCOLS_1D && ii~=size(this_data.c_all,1)
            figure;set(gcf,'Position',FIG_POS_1D)
        end
        
        clear preresp_idx thisrow thiscol;
    end
    
    
    % let's also try it the imagesc way...
    figure;
    for ii = 1:size(this_data.c_all,1)
        
        thisrow = 1;
        thiscol = this_data.r_all(ii);
        
        thiscol = mod(thiscol,NCOLS_1D);
        if thiscol==0
            thiscol = NCOLS_1D;
        end
        
        if thiscol==1 && this_data.r_all(ii)>NCOLS_1D && this_data.t_all(ii)==1
            % if this is the first trial of a new run on a new figure
            figure;set(gcf,'Position',FIG_POS_1D);
        end

        
        subplot(1,NCOLS_1D,(thisrow-1)*NCOLS_1D+thiscol); hold on;
        
        preresp_idx = ismember(this_data.s_all.XDAT{ii},[1 2]);

        
        imagesc(0.001*(1:sum(preresp_idx)),this_data.t_all(ii),sqrt( this_data.s_all.X{ii}(preresp_idx).^2 + this_data.s_all.Y{ii}(preresp_idx).^2 ).');
        
        if ~isempty(this_data.s_all.excl_trial{ii}) % if we do exclude
            % draw an x to the right of the trial
            plot(4.25,this_data.t_all(ii),'rx','LineWidth',1.5)
        end
        
        set(gca,'CLim',[0 1]*this_data.excl_criteria.delay_fix_thresh);

        axis ij off tight;
    end
    
    
    clear this_data;
    
end