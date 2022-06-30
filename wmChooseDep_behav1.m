% wmChooseDep_behav1.m
%
% behavioral stimulus presentation script for "choose the best" spatial WM
% MGS experiment
%
% participant sees 1 or 2 items (different colors, randomly chosen from a
% set of equiluminant colors), then remembers their positions precisely
% over delay (3.5 s). at end of delay, a response cue either tells
% participants to respond about a cued item (R1, R2_cued), like previous
% results (Sprague et al, 2014), or allows participants to choose an item
% (R2_choose). response is memory-guided saccade.
%
% TCS 6/24/2017
%
% TODO: stop recording w/ ESC
% TODO: fixation as dot within circle
%
% YML 6/1/2022
% TODO: pseudo randomize trial sequence for serial dependence


function wmChooseDep_behav1(subj,run)

try
p.expt_name = 'wmChooseDep_behav1';

p.do_et = 1;

p.subj = subj;
p.run = run;

% if ~exist('./data','dir')
%     mkdir('./data');
% end

p.filename = sprintf('../../data/wmChoose_data/%s_r%02.f_%s_%s.mat',p.subj,p.run,p.expt_name,datestr(now,30));
if p.do_et == 1
    p.eyedatafile = sprintf('%s_r%02.f',p.subj(1:min(length(p.subj),3)),p.run);
end

p.rng_seed = cputime*1000;
rng(p.rng_seed);


% ------ size of relevant stim features, etc ------ %
p.wm_ecc = 12;     % deg [in behavioral room, max ecc of circle 15 deg]
p.cue_size = 0.55; % deg
p.wm_size = 0.65;  % deg, size of WM dots
p.sep_ang = 30;    % deg polar angle, maximum stim separation distance
p.aperture_size = 15; % [or max ecc?]

p.fix_size_in  = 0.075; % radius, deg
p.fix_size_out = 0.30; % radius, deg
p.fix_pen = 1.5;
% now we draw out (frameoval), cue (if end of trial), in (filloval) for fix

% ------- keyboard stuff --------------------------- %
KbName('UnifyKeyNames');
if ismac == 1
    p.esc_key = KbName('escape'); % press this key to abort
else
    p.esc_key = KbName('escape');
end
%p.start_key = [KbName('5%') KbName('5')];  % should be lower-case t at prisma? (or %5, which is top-row 5, or 5, which is numpad 5)
p.space = KbName('space');


% ------ color of relevant stim features ----------- %
p.bg_color  = 20*[1 1 1];
p.fix_color = 75*[1 1 1];%[150 150 150];        % during trial/delay/etc
%p.wm_colors = [237 28 36;
%                 27 117 188;
%                 0 166 81;
%                 247 148 29;
%                 255 242 0;
%                 158 31 99;
%                 197 121 161;
%                 139 94 60];  % set of WM colors, equiluminant [TODO]

p.wm_colors = [200 0   0;         % red
                 0 0 255          % blue
                 180 0 180;       % purple
                 130 130 0];      % yellow
                 

p.choose_color = 130*[1 1 1];%[255 255 255]; % when subj should choose, color of fix



% ------ conditions ------ %
p.r_cond = [1 2 3]; % 1: R1, 2: R2_cued, 3: R2_choose
p.repetitions = 10; 
p.ntrials = length(p.r_cond)*p.repetitions;

p.conditions = nan(p.ntrials,1);
cnt = 1;
for cc = 1:length(p.r_cond)
    for rr = 1:p.repetitions
        p.conditions(cnt) = p.r_cond(cc);
        cnt = cnt+1;
    end
end

p.rnd_idx = randperm(p.ntrials);
p.conditions = p.conditions(p.rnd_idx,:);


% ------ timing of trial events --------- %
p.targ_dur = 0.5;
p.delay_dur = 3.5;
p.cue_dur = 1.5; % a bit longer than usual
p.feedback_dur = 0.8; 
p.iti_range = [2 4]; % randomly choose between those

p.itis = linspace(p.iti_range(1),p.iti_range(2),p.ntrials);
p.itis = p.itis(randperm(p.ntrials));



% ------- things to save -------- %
p.targ_coords = cell(2,1);
p.targ_coords{1} = nan(p.ntrials,2);
p.targ_coords{2} = nan(p.ntrials,2);

p.targ_colors = cell(2,1);
p.targ_colors{1} = nan(p.ntrials,3);
p.targ_colors{2} = nan(p.ntrials,3);

p.targ_angs = nan(p.ntrials,2);


% ------- Screen setup, optics --------- %

p.screen_height = 30; % cm, in the experiment room
p.viewing_distance = 56; % cm, in the experiment room (inside lab)

% open a screen, to get the resolution
[w, p.scr_rect] = Screen('OpenWindow',max(Screen('Screens')),[0 0 0]); HideCursor;
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);

p.ifi = Screen('GetFlipInterval',w);

p.center = p.scr_rect([3 4])/2;

p.ppd = p.scr_rect(4)/(2*atan2d(p.screen_height/2,p.viewing_distance));


p.aperture_rect = CenterRectOnPoint([0 0 2 2]*p.ppd*p.aperture_size,p.center(1),p.center(2));
p.fix_rect_out  = CenterRectOnPoint([0 0 2 2] * p.ppd  * p.fix_size_out,p.center(1),p.center(2));
p.fix_rect_in   = CenterRectOnPoint([0 0 2 2] * p.ppd  * p.fix_size_in, p.center(1),p.center(2));




% --------- eyetracking ----------- %
if p.do_et == 1
    
    el=EyelinkInitDefaults(w);
    
    el.backgroundcolour=p.bg_color(1);  % TODO: fix this?
    el.calibrationtargetcolour=p.fix_color(1);
    
    el.msgfontcolour=p.fix_color(1);
    p.foregroundcolour=p.fix_color(1);
    
    EyelinkUpdateDefaults(el);

    
    Eyelink('Initialize','PsychEyelinkDispatchCallback') % initialises the eyetracker
   
    Eyelink('command','calibration_type=HV13'); % updating number of callibration dots
    s=Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA');% (,GAZERES,HREF,PUPIL,STATUS,INPUT');
    s=Eyelink('command', 'sample_rate = 1000');
    s=Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
    s=Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
    
    
    
    % make sure that we get gaze data from the Eyelink
    % if a widescreen display:
    if p.scr_rect(3)==2560 && p.scr_rect(4)==1440
        %         Eyelink('command', 'generate_default_targets = NO');
        %
        %         % 13 samples
        %         Eyelink('command','calibration_samples = 13');
        %         Eyelink('command','validation_samples = 13');
        %
        %         % set up two random orders:
        %         calib_str = sprintf('%d,', randperm(13)-1);
        %         valid_str = sprintf('%d,', randperm(13)-1);
        %
        %         Eyelink('command',sprintf('calibration_sequence = %s',calib_str(1:end-1)));
        %
        %
        %         Eyelink('command',sprintf('validation_sequence = %s',valid_str(1:end-1)));
        
        Eyelink('command', 'calibration_area_proportion 0.59 0.83');
        %Eyelink('command', 'calibration_area_proportion 0.59 0.83');
        
    end
    
    
    %------ calibrate the eye tracker --------
    EyelinkDoTrackerSetup(el);
    if s~=0
        error('link_sample_data error, status: ',s)
    end
    Eyelink('openfile',p.eyedatafile);

end


% START OF EXPERIMENT

draw_aperture = @() Screen('FillOval',w,p.bg_color,p.aperture_rect);

Screen('FillRect',w,[0 0 0]);
draw_aperture();

% TODO: instructions/greeting
% instructions
txt = 'Remember dot position precisely';
Screen('TextSize', w, 30);
DrawFormattedText(w,txt,'center',p.center(2)-4*p.ppd,p.fix_color);


%Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
%Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2);
Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2);

Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2);
%Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
Screen('Flip',w);


% check for esc, space.... 

resp = 0;
while resp == 0
    [resp, ~] = checkForResp(p.space, p.esc_key);
    if resp == -1
        Screen('CloseAll'); ShowCursor;
        Eyelink('ShutDown'); 
        return;
    end
end
clear resp;
p.start_expt = GetSecs;


% blank screen
Screen('FillRect',w,[0 0 0]);
draw_aperture();
Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
%Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);

Screen('Flip',w);

if p.do_et == 1
    Eyelink('Message','xDAT %i', 10);
    Eyelink('StartRecording'); % make 1 big edf file (save time)
end


WaitSecs(1.5); % wait a bit before first trial


for tt = 1:p.ntrials
    
    
    % this trial's position(s)
    this_ang = nan(1,2);
    this_ang(1) = 360*rand(1);
    
    
    p.targ_coords{1}(tt,:) = p.wm_ecc * [cosd(this_ang(1)) sind(this_ang(1))];
    
    if p.conditions(tt,1) ~= 1
        this_ang(2) = this_ang(1)+p.sep_ang + (360-2*p.sep_ang)*rand(1); % compute second position if necessary
        p.targ_coords{2}(tt,:) = p.wm_ecc * [cosd(this_ang(2)) sind(this_ang(2))];
    end
    
    p.targ_angs(tt,:) = this_ang; % save these for convenience
    
    
    % this trial's colors
    tmp_color_idx = randperm(size(p.wm_colors,1));
    
    p.targ_colors{1}(tt,:) = p.wm_colors(tmp_color_idx(1),:);
    
    if p.conditions(tt,1) ~= 1
        p.targ_colors{2}(tt,:) = p.wm_colors(tmp_color_idx(2),:);
    end
    
    
    
    % targets (XDAT 1) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    trial_start = GetSecs;
    
    
    if p.do_et == 1
        %Eyelink('Message','TarX %s', num2str(p.targ_coords{1}(tt,1)));
        %Eyelink('Message','TarY %s', num2str(p.targ_coords{1}(tt,2)));
        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));
        
        
        Eyelink('Message','xDAT %i',1);
        
        Eyelink('command', 'record_status_message "TRIAL %d of %d"', tt, p.ntrials);
        
    end
    
    
    
    while GetSecs < trial_start + p.targ_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        % target 1
        Screen('DrawDots',w,p.ppd*[1;-1].*p.targ_coords{1}(tt,:).', p.wm_size*p.ppd, p.targ_colors{1}(tt,:), p.center, 2);
        
        if p.conditions(tt,1)~=1
            
            % if necessary, target 2
            Screen('DrawDots',w,p.ppd*[1;-1].*p.targ_coords{2}(tt,:).', p.wm_size*p.ppd, p.targ_colors{2}(tt,:), p.center, 2);
        end
        
        % fixation
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
        
    end
    
    % delay (XDAT 2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    if p.do_et == 1
        Eyelink('Message','xDAT %i',2);    
    end
    
    
    while GetSecs < trial_start + p.targ_dur + p.delay_dur
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
%         if p.conditions(tt,1) ~= 3
%             Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.choose_color,p.center,2);
%         else
%             Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.targ_colors{1}(tt,:),p.center,2);
%         end
        
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
        
    end
    
    % go cue (XDAT 3) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1
        Eyelink('Message','xDAT %i',3);
    end
    
    while GetSecs < trial_start + p.targ_dur + p.delay_dur + p.cue_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        if p.conditions(tt,1) == 3
            Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.choose_color,p.center,2);
        else
            Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.targ_colors{1}(tt,:),p.center,2);
        end
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
    end
    
    % feedback (XDAT 4, tarx, tary) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(p.targ_coords{1}(tt,1)));
        Eyelink('Message','TarY %s', num2str(p.targ_coords{1}(tt,2)));
        % NOTE: incorrect on 1/3 trials
        Eyelink('Message','xDAT %i',4);
        
    end
    
    while GetSecs < trial_start + p.targ_dur + p.delay_dur + p.cue_dur + p.feedback_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        % target 1
        Screen('DrawDots',w,p.ppd*[1;-1].*p.targ_coords{1}(tt,:).', p.wm_size*p.ppd, p.targ_colors{1}(tt,:), p.center, 2);
        
        if p.conditions(tt,1)~=1
            
            % if necessary, target 2
            Screen('DrawDots',w,p.ppd*[1;-1].*p.targ_coords{2}(tt,:).', p.wm_size*p.ppd, p.targ_colors{2}(tt,:), p.center, 2);
        end
        
        % fixation
%        Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

        
         %       Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        if p.conditions(tt,1) == 3
            Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.choose_color,p.center,2);
        else
            Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.targ_colors{1}(tt,:),p.center,2);
        end
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

        
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
 
    end
    
    % ITI (XDAT 5) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));
        
        Eyelink('Message','xDAT %i',5);
    end

    while GetSecs < trial_start + p.targ_dur + p.delay_dur + p.cue_dur + p.feedback_dur + p.itis(tt)
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');

            save(p.filename,'p');
            return;
        end
        
    end
    
    
    % save [note: in scanner, do this at beginning of ITI after first flip]
    save(p.filename,'p');
    
end

% END OF EXPERIMENT - TEXT
if p.do_et == 1
    Eyelink('Message','xDAT %i',11);
end

Screen('FillRect',w,[0 0 0]);
draw_aperture();
txt = sprintf('End of run %i',p.run);
DrawFormattedText(w,txt,'center',p.center(2)-4*p.ppd,p.fix_color);
%Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);

Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

Screen('Flip',w);

resp = 0;
while resp == 0
    [resp, timeStamp] = checkForResp(p.space, p.esc_key);
end
clear resp;

if p.do_et == 1
    Eyelink('StopRecording');
    Eyelink('ReceiveFile',[p.eyedatafile '.edf'],[p.eyedatafile '.edf']);
    
    p.eyedatafile_renamed = [p.filename(1:(end-3)) 'edf'];
    movefile([p.eyedatafile '.edf'],p.eyedatafile_renamed);
    
    Eyelink('ShutDown');
end

Screen('CloseAll');
ShowCursor;
catch
    Screen('CloseAll');
    ShowCursor;
end


return