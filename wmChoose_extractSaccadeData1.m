% wmChoose_extractSaccadeData1.m (from MGSMap_extractSaccadeData1.m)
%
% extracts simple versions of initial, final saccades, saves them for use





function wmChoose_extractSaccadeData1(subj)

if nargin < 1 || isempty(subj)
    %subj = {'KD','CC','EK','MR','AB'};
    subj = {'aa','ab1','ac1','ac2','ae','af','ag'}; %aa1

end

% if nargin < 2 || isempty(sess)
%     sess = {{'MGSMap1','MGSMap2'},{'MGSMap1','MGSMap2'},{'MGSMap1'},{'MGSMap1'},{'MGSMap1'}};
% end

root = '/Volumes/data/wmChoose';

for ss = 1:length(subj)
    %for thissess = 1:length(sess{ss})
        
        % load already-concatenated behavioral data
        behav_fn = sprintf('%s/data/%s_wmChoose_behav.mat',root,subj{ss});
        thisbehav = load(behav_fn);
        
        % load preproc'd eye data
        iEye_fs = sprintf('%s/preproc_iEye/%s_wmChoose_behav1_r*_preproc.mat',root,subj{ss});
        iEye_fn = dir(iEye_fs);
        
        
        % for all iEye files...
        for ii = 1:length(iEye_fn)
            
            fprintf('loading %s/%s\n',iEye_fn(ii).folder,iEye_fn(ii).name);
            this_et = load(sprintf('%s/%s',iEye_fn(ii).folder,iEye_fn(ii).name));
        
            this_cfg = this_et.ii_cfg;
            this_data = this_et.ii_data;
            
            if ii == 1
                
%                 nblank = length(iEye_fn)*this_cfg.numtrials;
%                 
%                 i_sacc = nan(nblank,2);
%                 f_sacc = nan(nblank,2);
%                 
%                 i_sacc_err = nan(nblank,1);
%                 f_sacc_err = nan(nblank,1);
%                 
%                 i_sacc_rt = nan(nblank,1);

                s_all = thisbehav.s_all;
                
                % also want to store the raw coordinates
                s_all.i_sacc_raw = nan(size(s_all.i_sacc));
                s_all.f_sacc_raw = nan(size(s_all.f_sacc));
                
                s_all.i_sacc = nan(size(s_all.i_sacc));
                s_all.f_sacc = nan(size(s_all.f_sacc));
                
                
                s_all.i_sacc_err = nan(size(s_all.i_sacc,1),1);
                s_all.f_sacc_err = nan(size(s_all.f_sacc,1),1);
                
                s_all.ti_all = nan(size(s_all.i_sacc,1),size(this_cfg.trialinfo,2));
                
                
            end
            
            thisidx = (1:this_cfg.numtrials)+(ii-1)*this_cfg.numtrials;
            
            % % % %
            % trialinfo contains:
            % - condition - (1) R1, (2) R2, (3) R2-chooose
            % - coords_1 (2 cols)
            % - coords_2 (2 cols, nan if R1)
            % - 'chosen' item (1 or 2)
            % - distance from each of the two targets (2 cols)
            
            % XDAT: 1 = targ, 2 = delay, 3 = go, 4 = feedback, 5
            
            
            
            % extract initial saccade (first saccade after go cue)
        
            [this_data,this_cfg] = ii_selectfixationsbytrial(this_data,this_cfg,'XDAT',3,'first');
            
            % in case there aren't n_trials fixations, have to look up
            % which trial goes where
            for tt = 1:size(this_cfg.cursel,1)
                s_all.i_sacc_raw(thisidx(this_cfg.trialvec(this_cfg.cursel(tt,1))),:) = [this_data.X_fix(this_cfg.cursel(tt,1)) this_data.Y_fix(this_cfg.cursel(tt,1))];
            end
        
            % extract final saccade (last saccade before feedback stim)
        
            [this_data,this_cfg] = ii_selectfixationsbytrial(this_data,this_cfg,'XDAT',3,'last');
            %s_map.f_sacc(thisidx,:) = [this_data.X_fix(this_cfg.cursel(:,1)) this_data.Y_fix(this_cfg.cursel(:,1))];
            % in case there aren't n_trials fixations, have to look up
            % which trial goes where
            for tt = 1:size(this_cfg.cursel,1)
                s_all.f_sacc_raw(thisidx(this_cfg.trialvec(this_cfg.cursel(tt,1))),:) = [this_data.X_fix(this_cfg.cursel(tt,1)) this_data.Y_fix(this_cfg.cursel(tt,1))];
            end
            
            % compute error for initial; final saccades given which
            % stimulus is 'chosen' [[maybe do both??]]
            % error for initial
            for tt = 1:length(thisidx)
                
                this_coord = this_cfg.trialinfo(tt,[1 2] + 1 + (this_cfg.trialinfo(tt,6)-1)*2);
                
                s_all.i_sacc_err(thisidx(tt)) = sqrt(sum((s_all.i_sacc_raw(thisidx(tt),:)-this_coord).^2,2));
                s_all.f_sacc_err(thisidx(tt)) = sqrt(sum((s_all.f_sacc_raw(thisidx(tt),:)-this_coord).^2,2));
            end

            s_all.ti_all(thisidx,:) = this_cfg.trialinfo; % save this too, which includes which item was chosen
        
            
            clear thisidx this_data this_cfg this_et;
        
        end

        % need to loop over trials because we have to read out the chosen
        % coord
        
        % if NAN, do nothing.... (we'll have to cull those trials
        % elsewhere?)
        
        for tt = 1:size(s_all.i_sacc_raw,1)
            
            if ~isnan(s_all.ti_all(tt,6))
                % rotate all so that i_sacc, f_sacc are at known position...
                [tmpth,tmpr] = cart2pol(s_all.i_sacc_raw(tt,1),s_all.i_sacc_raw(tt,2));
                tmpth = tmpth - deg2rad(thisbehav.c_all(tt,1+s_all.ti_all(tt,6)));
                [s_all.i_sacc(tt,1),s_all.i_sacc(tt,2)] = pol2cart(tmpth,tmpr);
                clear tmpth tmpr;
                
                [tmpth,tmpr] = cart2pol(s_all.f_sacc_raw(tt,1),s_all.f_sacc_raw(tt,2));
                tmpth = tmpth - deg2rad(thisbehav.c_all(tt,1+s_all.ti_all(tt,6)));
                [s_all.f_sacc(tt,1),s_all.f_sacc(tt,2)] = pol2cart(tmpth,tmpr);
                clear tmpth tmpr;
            end
        end
        
        % append s_map to behavioral data file
        save(behav_fn,'s_all','-append')
        
        clear s_map thisbehav behav_fn iEye_fn;
        
    %end
end





return