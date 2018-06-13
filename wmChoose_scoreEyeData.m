% wmChoose_scoreEyeData.m
%
% preprocesses, scores, and QC's data using new autoscore in iEye_ts
% - replaces wmchoose_runPreProc.m, wmChoose_extractSaccadeData.m, and
%   wmChoose_plotSummary.m

close all;
root = '/Volumes/data/wmChoose';
ifg_fn = '~/Documents/MATLAB/toolboxes_dev/iEye_ts/examples/p_1000hz.ifg';

subj = {'aa1','aa2','ab1','ab2','ac1','ac2','ae','af','ag','ai'}; %aa1
%subj = {'ah'};

runs_with_err = {};
errs = {};

fn_prefix = 'wmChoose_behav1'; % OR _ds_preCue

QCdir = fullfile(root,'scoring_QC');

WHICH_EXCL = [11 13 20 21 22]; % based on most recent version of plotResults

% set up iEye params
ii_params = ii_loadparams;
ii_params.trial_end_value = 5;
ii_params.drift_epoch = [1 2];
ii_params.calibrate_epoch = 4;
ii_params.calibrate_select_mode = 'last';
ii_params.calibrate_window = 300; % originally was 300
ii_params.blink_window = [200 200];
ii_params.plot_epoch = [3 4];
ii_params.calibrate_limits = [2.5]; % error during feedback shouldn't exceed this

ii_params.ppd = 34.1445; % behavioral room screen


for ss = 1:length(subj)
    
    fns = sprintf('%s/data/%s_r*_%s_*.edf',root,subj{ss},fn_prefix);
    thisf = dir(fns);
    clear fns;
    
    % remove files that are "r00"
    ii_trial = cell(length(thisf),1);
    block_num = nan(length(thisf),1);
    
    for ff = 1:length(thisf)
        
        fprintf('Preprocessing %s\n',thisf(ff).name);
        
        this_edf = sprintf('%s/data/%s',root,thisf(ff).name);
        
        
        % look for matching mat file
        %fns = sprintf('%s/data/%s_%s*block%02.f_*.mat',root,subj{ss},fn_prefix,block_num);
        %matf = dir(fns);
        matf = sprintf('%smat',this_edf(1:end-3));
        thisbehav = load(matf);

        % for convenience...
        thisbehav = thisbehav.p;
        
        % custom for each expt
        block_num(ff) = str2double(matf(strfind(matf,'_r')+[2 3]));
        
        % turn targ_coords variable into cell of cells
        % TODO: allow each cell to have numel==n_selections * 2
        coords = cell(thisbehav.ntrials,1);
        for tt = 1:thisbehav.ntrials
            coords{tt} = {thisbehav.targ_coords{1}(tt,:), thisbehav.targ_coords{2}(tt,:)};
        end
        
        
        % set up trialinfo
        trial_info = horzcat(thisbehav.conditions,thisbehav.targ_coords{:});


        
        preproc_fn = sprintf('%s/preproc_iEye/%s_%s_r%02.f_preproc.mat',root,subj{ss},fn_prefix,block_num(ff));
        
       
        [ii_data,ii_cfg,ii_sacc] = wmChoose_preproc1(this_edf,ifg_fn,preproc_fn,coords,trial_info,ii_params);
        
        this_targ_coords = nan(size(ii_cfg.trialinfo,1),2);
        this_targ_coords(ismember(ii_cfg.trialinfo(:,1),[1 2]),:) = ii_cfg.trialinfo(ismember(ii_cfg.trialinfo(:,1),[1 2]),[2 3]);
        
        thisidx = find(ii_cfg.trialinfo(:,1)==3);
        for tt = 1:length(thisidx)
            this_targ_coords(thisidx(tt),:) = ii_cfg.trialinfo(thisidx(tt),ii_cfg.trialinfo(thisidx(tt),6)*2+[0 1]);
        end
        
        [ii_trial{ff}, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,this_targ_coords,3,ii_params.drift_epoch);
        
        
         
        close all;
        close all hidden;
        clear preproc_fn trial_info cond thisbehav matf fns this_edf thisbehav this_targ_coords;
    end
    
    % save for this subj a wmChoose_scored.mat containing concatenated
    % ii_trial (ii_sess) and ii_cfg
    ii_sess = ii_combineruns(ii_trial,block_num);
    
    save(sprintf('%s/data/%s_wmChoose_scored.mat',root,subj{ss}),'ii_sess','WHICH_EXCL');
    
    % plot QC files and save them [TODO]
    
    % exclusion report
    fh_excl = ii_plotQC_exclusions(ii_sess,ii_cfg,WHICH_EXCL,0);
    for ff = 1:length(fh_excl)-1
        saveas(fh_excl(ff),sprintf('%s/%s_excl_dot_%02.f.png',QCdir,subj{ss},ff));
    end
    saveas(fh_excl(end),sprintf('%s/%s_excl_all.png',QCdir,subj{ss}),0);
    close(fh_excl);
    
    % all trials
    fh_trials = ii_plotQC_alltrials(ii_sess,ii_cfg,WHICH_EXCL,0);
    for ff = 1:length(fh_trials)
        saveas(fh_trials(ff),sprintf('%s/%s_trials_r%02.f.png',QCdir,subj{ss},block_num(ff)));
    end
    close(fh_trials);
    
    close all hidden;
    clear ii_sess;
end