% ssPri_concatBehav1_map.m
%
% combines all the behavioral data files into one for easy loading
%
% TCS, adapted from ssPri_concatBehav1.m, used for mapping-only scans,
% should they exist (e.g., CCm)
% (note - not actually necessary, if no main-task files found, then nothing
% is done w/ them)
%
% TCS updating 3/5/2017 to allow for missing eyetracker data on some runs
%
% TCS updated 4/19/2017 to conform to new data organization in
% ssPri_scanner_afni; subj/sess breakdown
%
%
% c_all: 1) cue conditions (R1, R2-cue, R2-choose)
%        2) remembered angle 1 (cued, for R1/R2-cue)
%        3) remembered angle 2 (uncued for R1/R2-cue)

root_behav  = 'Z:/projects/wmChoose/data';%CC_MGSMap25mm_MB4_behav'; %HACK
root_target = 'Z:/projects/wmChoose/data'; % where to save things

subj = {'sub001','sub002','sub003','sub004','sub005','sub006','sub007','sub008','sub009','sub010','sub011','sub012','sub013','sub014','sub016','sub017','sub018','sub019','sub020','sub021'};
%subj = {'ah','ai'};
%sess = {{'MGSMap2'}};%,{'Map1','Map2'},{'Map1','Map2'},{'Map1','Map2'}}; % need to turn this into KDm1, KDm2 for behavioral data extraction

%n_subRuns = 2; % could do this automagically...


for ss = 1:length(subj)
    
    %for sess_idx = 1:length(sess{ss})
        
        this_subjID = sprintf('%s',subj{ss});
        
%         if strfind(sess{ss}{sess_idx},'Map')
%             this_subjID = [subj{ss} strrep(sess{ss}{sess_idx},'Map','m')];
%         else
%             this_subjID = [subj{ss} sess{ss}{sess_idx}]; % should never use this...
%         end
        % start w/ mapping data
        
        % look for all behavioral files (not saccades yet) (like: ag_r09_wmChoose_behav1_20170822T141641.mat

        fnm_b = sprintf('%s/%s_r*_wmChoose_behav1_*.mat',root_behav,this_subjID); % because of datestr(now,30), there's a T in the behav files, but not in the saccade files
        tmp_fm_b = dir(fnm_b);
        
        % get rid of "_iEye" files if they exist
        tmp_iEye = strfind({tmp_fm_b.name},'iEye');
        myfm_b = tmp_fm_b(cellfun(@isempty,tmp_iEye));
        
        clear tmp_iEye tmp_fm_b;
        
        % myfm_b
        
        % and look for saccade files (loop through behav files, look for each
        % saccade file, otherwise store an empty fname)
%         for ii = 1:length(myfm_b)
%             fnm_s = sprintf('%s/%s_ssPri_ds_preCue_scanMapping_run%02.f_*saccadeData.mat',root_behav,this_subjID,ii); % because of datestr(now,30), there's a T in the behav files, but not in the saccade files
%             this_f = dir(fnm_s);
%             if ~isempty(this_f)
%                 myfm_s(ii) = dir(fnm_s);
%             else myfm_s(ii).name = [];
%             end
%             clear this_f;
%         end
%         
        % initialize the variables we'll fill up
        %t_map = [];   % timing
        c_all = [];   % conditions/trial labels
        p_all.radMean = []; % stimulus parameters
        p_all.polarAngleOffset = [];
        r_all = [];  % run #
        t_all = []; % trial #
        
        coords_all{1} = [];
        coords_all{2} = [];
        colors_all{1} = [];
        colors_all{2} = [];
        
        % (these should be same length - should probably check for that here)
        for ff = 1:length(myfm_b)
            
            fnb = sprintf('%s/%s',root_behav,myfm_b(ff).name);
            fprintf('loading %s...\n',fnb);
            bdata = load(fnb);
            
            rnum = str2double(myfm_b(ff).name(strfind(myfm_b(ff).name,'_r')+[2 3]));
            
            r_all = [r_all; rnum*ones(bdata.p.ntrials,1)];
            t_all = [t_all; (1:bdata.p.ntrials).'];
            %t_map = [t_map; bdata.tasktiming.segStart - bdata.tasktiming.expStartTime];
            %t_map = [t_map; [bdata.p.delay_start bdata.p.resp_start]-bdata.p.expt_start];
            c_all = [c_all; bdata.p.conditions bdata.p.targ_angs];
            p_all.radMean = [bdata.p.wm_ecc];
            
            coords_all{1} = [coords_all{1}; bdata.p.targ_coords{1}];
            coords_all{2} = [coords_all{2}; bdata.p.targ_coords{2}];

            colors_all{1} = [colors_all{1}; bdata.p.targ_colors{1}];
            colors_all{2} = [colors_all{2}; bdata.p.targ_colors{2}];

            
            %p_all.polarAngleOffset = [bdata.p.offset];
            
            
%             if ~isempty(myfm_s(ff).name)
%                 fns = sprintf('%s%s/%s',root_behav,this_subjID,myfm_s(ff).name);
%                 fprintf('loading %s...\n',fns);
%                 sdata = load(fns);
%                 
%                 % saccade data
%                 s_all.i_sacc = [s_all.i_sacc;sdata.i_sacc_coord];
%                 s_all.f_sacc = [s_all.f_sacc; sdata.f_sacc_coord];
%                 s_all.i_sacc_err = [s_all.i_sacc_err; sdata.i_sacc_err];
%                 s_all.f_sacc_err = [s_all.f_sacc_err; sdata.f_sacc_err];
%                 
%                 % RTs
%                 thissrt = nan(size(sdata.sacc_onset,1),1);
%                 thistrials = ~isnan(sdata.i_sacc);
%                 thissrt(thistrials) = sdata.sacc_onset(sub2ind(size(sdata.sacc_onset),find(thistrials),sdata.i_sacc(thistrials.')));
%                 s_all.srt = [s_all.srt; 2*thissrt];% sdata.sacc_onset(sub2ind(size(sdata.sacc_onset),1:size(sdata.sacc_onset,1),sdata.i_sacc.'))  ];
%                 % NOTE: because of 500 Hz eyetrakcing in scanner, using 2*
%                 % rt above
%                 
%                 clear sdata thissrt thistrials;
%             %else
            %end
            
            
            % WE WANT TO SAVE:
            % - timing (bdata.tasktiming.segStart -
            %   bdata.tasktiming.expStartTime
            % - conditions - [bdata.task.polarAng.' bdata.task.targQuads]
            % - stimulus parameters (mapping) - bdata.stimulus.radMean
            % - saccade endpoint(s)
            % - saccade RTs
            
            
            
            clear bdata;
        end
        
        
        
        
        % save everything
        fn2s = sprintf('%s/%s_wmChoose_behav.mat',root_target,subj{ss});
        fprintf('saving to %s...\n',fn2s);
        save(fn2s,'c_all','p_all','r_all','t_all','coords_all','colors_all');
        
        % clear those things..
    %end
end