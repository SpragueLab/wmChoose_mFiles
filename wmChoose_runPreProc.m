
close all;
root = '/Volumes/data/wmChoose';
ifg_fn = '/Volumes/tommy/Documents/MATLAB/toolboxes_dev/iEye_ts/examples/p_1000hz.ifg';

%subj = {'aa1','aa2','ab1','ab2','ac1','ac2','ae','af','ag'}; %aa1
subj = {'ae'};

runs_with_err = {};
errs = {};

fn_prefix = 'wmChoose_behav1'; % OR _ds_preCue

% set up iEye params
ii_params = ii_loadparams;
ii_params.trial_end_value = 5;
ii_params.drift_epoch = [1 2];
ii_params.calibrate_epoch = 4;
ii_params.calibrate_select_mode = 'nearest';
ii_params.blink_window = [200 200];
ii_params.plot_epoch = [3 4];
ii_params.calibrate_limits = [0.75 1.333]; % original ecc b/w 9 and 16...

ii_params.ppd = 34.1445; % behavioral room screen


for ss = 1:length(subj)
    
    fns = sprintf('%s/data/%s_r*_%s_*.edf',root,subj{ss},fn_prefix);
    thisf = dir(fns);
    clear fns;
    
    % remove files that are "r00"
    
    
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
        block_num = str2double(matf(strfind(matf,'_r')+[2 3]));
        
        % turn targ_coords variable into cell of cells
        % TODO: allow each cell to have numel==n_selections * 2
        coords = cell(thisbehav.ntrials,1);
        for tt = 1:thisbehav.ntrials
            coords{tt} = {thisbehav.targ_coords{1}(tt,:), thisbehav.targ_coords{2}(tt,:)};
        end
        
        
        % set up trialinfo
        trial_info = horzcat(thisbehav.conditions,thisbehav.targ_coords{:});


        
        preproc_fn = sprintf('%s/preproc_iEye/%s_%s_r%02.f_preproc.mat',root,subj{ss},fn_prefix,block_num);
        
       % try
            wmChoose_preproc1(this_edf,ifg_fn,preproc_fn,coords,trial_info,ii_params);
%         catch
%             runs_with_err{end+1} = this_edf;
%             errs{end+1} = lasterror;
%         end
        
        close all;
        close all hidden;
        clear preproc_fn trial_info cond thisbehav matf fns block_num this_edf thisbehav;
    end
    
end