function wmPri_scoreEyeData(subj,sess)
% preprocesses eye position data (EDF) files from wmPri scanning experiment
%
% TCS 8/17/2018

% feedback = 6
% pre-targets = 1
%
% trial start/end = 1/7
%

close all;
root = wmPri_loadRoot;
ifg_fn = 'Y:/matlab/toolboxes/iEye_ts/examples/p_500hz.ifg';

WHICH_EXCL = [11 13 20 21 22]; % everything except calibration errors (for now)

if nargin < 1
    %subj = {'JK','YK','MH'};
    subj = {'MH'};
    
end

if nargin < 2
    %sess = {{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1'}};
    
    
end

QCdir = 'Z:/projects/nyu/wmPri/wmPri_iEye_score_QC';

fn_prefix = 'wmPri_scanner_sp'; % OR _ds_preCue

% set up iEye params
ii_params = ii_loadparams;
ii_params.trial_end_value = 7;
ii_params.drift_epoch = [1 2 3 4];
ii_params.calibrate_epoch = 6;
ii_params.calibrate_mode = 'run';
ii_params.calibrate_select_mode = 'nearest';
ii_params.blink_window = [200 200];
ii_params.blink_thresh = 3.5; %default: 1.5 (percentile)
ii_params.plot_epoch = [4 5 6];
ii_params.calibrate_limits = 2.5; % original ecc b/w 9 and 16...

ii_params.resolution = [1280 1024];

ii_params.ppd = 32.981; % for scanner, 1280 x 1024

for ss = 1:length(subj)
    
    for sessidx = 1:length(sess{ss})
        
        fns = sprintf('%s/rawbehav/%s_%s_behav/%s_%s_r*_%s_*.edf',root,subj{ss},sess{ss}{sessidx},subj{ss},sess{ss}{sessidx},fn_prefix);
        thisf = dir(fns);
        clear fns;
        
        run_labels = nan(length(thisf),1);
        ii_trial = cell(length(thisf),1);
        
        for ff = 1:length(thisf)
            
            fprintf('Preprocessing %s\n',thisf(ff).name);
            
            this_edf = sprintf('%s/rawbehav/%s_%s_behav/%s',root,subj{ss},sess{ss}{sessidx},thisf(ff).name);
            
            
            % look for matching mat file
            %fns = sprintf('%s/data/%s_%s*block%02.f_*.mat',root,subj{ss},fn_prefix,block_num);
            %matf = dir(fns);
            matf = sprintf('%smat',this_edf(1:end-3));
            thisbehav = load(matf);
            
            % for convenience...
            thisbehav = thisbehav.p;
            
            % custom for each expt
            block_num = str2double(matf(strfind(matf,'_r')+[2 3]));
            run_labels(ff) = block_num;
            
            % set up trialinfo
            % - to match w/ wmChoose, we'll store condition #, then xy for
            % item 1 (cued) and item 2 (uncued)
            
            coords = cell(thisbehav.ntrials,1);
            for tt = 1:thisbehav.ntrials
                coords{tt} = {thisbehav.targ_coords{1}(tt,:), thisbehav.targ_coords{2}(tt,:)};
            end
            
            
            % set up trialinfo
            trial_info = horzcat(thisbehav.conditions,thisbehav.targ_coords{:});
            
            
            
            preproc_fn = sprintf('%s/wmPri_iEye_preproc/%s_%s_r%02.f_preproc.mat',root,subj{ss},sess{ss}{sessidx},block_num);
            
            % preprocess raw data and extract saccades
            [ii_data,ii_cfg,ii_sacc] = ii_preproc(this_edf,ifg_fn,preproc_fn,ii_params,trial_info);
            
            % define resp epoch, fix epoch, etc. note that default behavior
            % of ii_scoreMGS will pick these up automatically
            ii_trial{ff} = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,[],ii_params.calibrate_epoch-1,ii_params.drift_epoch);
            
            close all;
            clear preproc_fn trial_info cond thisbehav matf fns block_num this_edf thisbehav;
        end
        
        ii_sess = ii_combineruns(ii_trial,run_labels);
        
        save(sprintf('%s/wmPri_behav/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx}),'ii_sess','WHICH_EXCL');

        
        % exclusion report
        fh_excl = ii_plotQC_exclusions(ii_sess,ii_cfg,WHICH_EXCL,0);
        for ff = 1:length(fh_excl)-1
            saveas(fh_excl(ff),sprintf('%s/%s_%s_excl_dot_%02.f.png',QCdir,subj{ss},sess{ss}{sessidx},ff));
        end
        saveas(fh_excl(end),sprintf('%s/%s_%s_excl_all.png',QCdir,subj{ss},sess{ss}{sessidx}),0);
        close(fh_excl);
        
        % all trials
        fh_trials = ii_plotQC_alltrials(ii_sess,ii_cfg,WHICH_EXCL,0);
        for ff = 1:length(fh_trials)
            set(fh_trials(ff),'Renderer','painters'); % hack to make sure saving as png works...
            saveas(fh_trials(ff),sprintf('%s/%s_%s_trials_r%02.f.png',QCdir,subj{ss},sess{ss}{sessidx},run_labels(ff)));
        end
        close(fh_trials);
        
        close all hidden;
        clear ii_sess;

        
    end
    
end


end