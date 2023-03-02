%% DoG
% c = sqrt(2)/exp(-0.5);
% a = 2;
% w = 1;
% 
% x = [-180:20:180];
% 
% y = x.*a.*w.*c.*exp(-(w.*x).^2);
% 
% figure; hold on
% plot(x,y)
% 
% % c = sqrt(2)/exp(-0.5);
% a = -6;
% w = 3;
% 
% % x = [-4:0.005:4];
% 
% y = x.*a.*w.*c.*exp(-(w.*x).^2);
% 
% % figure;
% plot(x,y)

%% Call in data

root = 'Z:/projects/wmChooseSD1';
% subj = {'sub001','sub002','sub003','sub020','sub022','sub024','sub025','sub026','sub027','sub028','sub032','sub034','sub035','sub036','sub037','sub039','sub043','sub044','sub045','sub046','sub047'};
subj = {'sub001','sub002','sub003','sub020','sub022','sub024','sub025','sub026','sub027','sub028','sub032','sub034','sub035','sub036','sub037','sub044','sub045','sub046','sub047'};
% subj = {'sub001'};
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
    %fn = sprintf('%s/data/%s_wmChoose_behav.mat',root,subj{ss});
    fn = sprintf('%s/data/%s_wmChooseSD1_behav.mat',root,subj{ss});
    fprintf('Loading trial information from %s\n',fn);
    this_data = load(fn);
    
    fn = sprintf('%s/data/%s_wmChooseSD1_scored.mat',root,subj{ss});
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


%% table of % included or excluded trials from each subj

percent_useTrial = nan(length(u_subj),1);
percent_excTrial = nan(length(u_subj),1);
num_useTrial = nan(length(u_subj),1);


for ss = 1:length(u_subj)

    thisidx = all_data.subj_all==ss & all_data.use_trial==1;
    thisidx2 = all_data.subj_all==ss & all_data.use_trial==0;
    thisidx3 = all_data.subj_all==ss;
    percent_useTrial(ss) = sum(thisidx)/sum(thisidx3);
    percent_excTrial(ss) = sum(thisidx2)/sum(thisidx3);
    num_useTrial(ss) = sum(thisidx);

end


%% Prep variables & param
c = sqrt(2)/exp(-0.5);

%convert all xy coordinate to angle
all_data.s_all.targAng360 = mod(atan2d(all_data.s_all.targ(:,2),all_data.s_all.targ(:,1)),360);
%all_data.s_all.f_sacc_RawAng = atan2d(all_data.s_all.f_sacc_raw(:,2),all_data.s_all.f_sacc_raw(:,1));
all_data.s_all.f_sacc_ang360 = mod(atan2d(all_data.s_all.f_sacc_raw(:,2),all_data.s_all.f_sacc_raw(:,1)),360); %all target is align to one point like 0 degree, i.e., reported loc - target


%compute target difference with current trial target - previous trial
%convert all xy coordinate to angle
all_data.s_all.targAng180 = atan2d(all_data.s_all.targ(:,2),all_data.s_all.targ(:,1));
all_data.s_all.f_sacc_RawAng = atan2d(all_data.s_all.f_sacc_raw(:,2),all_data.s_all.f_sacc_raw(:,1));
all_data.s_all.f_sacc_align180 = atan2d(all_data.s_all.f_sacc(:,2),all_data.s_all.f_sacc(:,1)); %all target is align to one point like 0 degree, i.e., reported loc - target

targ_diff = nan(length(all_data.t_all),1);
% targ_diff2 = nan(length(all_data.t_all) - length(all_data.r_all),1);

for ss = 1:length(u_subj)

    
    ru = unique(all_data.r_all(all_data.subj_all == ss));
    tu = unique(all_data.t_all(all_data.subj_all == ss));

    for rr = 1:length(ru)
        for tt = 2:length(tu)

            thisidx1 = all_data.subj_all==ss & all_data.r_all == ru(rr) & all_data.t_all == tt; %without ru or tu sometimes cause error with consecutive
            thisidx2 = all_data.subj_all==ss & all_data.r_all == ru(rr) & all_data.t_all == tt-1;

            targ_diff(thisidx1) = angdiffdeg(all_data.s_all.targAng180(thisidx1),all_data.s_all.targAng180(thisidx2)); %current trial ang - previous trial ang

        end
    end
end
all_data.use_trial(all_data.t_all == 1) = 0;

%% plot target loc against report loc

figure;
for ss = 1:length(u_subj)
    subplot(3,7,ss); hold on
    histogram(all_data.s_all.f_sacc_err, 45,'BinLimits',[0 5]);

    xlabel('Final Saccade Error (MAE;°)');
    title(u_subj{ss});

end



figure;
for ss = 1:length(u_subj)


        subplot(3,7,ss); hold on;

        thisidx = all_data.subj_all==ss & all_data.use_trial==1;
        thisidx2 = find(thisidx);
        
   
           scatter(all_data.s_all.targAng360(thisidx2),all_data.s_all.f_sacc_ang360(thisidx2),10,'r','filled', 'markerfacealpha',.5);
               

     xlim([0 360]);ylim([0 360]);
     set(gca,'TickDir','out','XTick',[0 180 360],'YTick',[0 180 360]); 
     xlabel('Target location (°)'); ylabel('Reported location (°)');
     axis square;
    
end

% %% plot error as function of targ location
% my_Bin = 0:30:360; % create bin limit for 0 to 360, make 1 bin for 30 degree interval therefore result in 12 bins in total
% %[BINS, EDGES] = discretize(X,N) where N is a scalar integer, divides the range of X into N uniform bins, and also returns the bin edges.
% location_Bin = discretize(all_data.s_all.targAng360,my_Bin);%put target location data into bin from 1-12, e.g., 1 is 0-30 degree and so on.
% 
% bin_colors = hsv(length(my_Bin));%give each bin a distinct color
% 
% figure;
% for ss = 1:length(u_subj) % go through all the subjects
% 
% 
%     subplot(3,7,ss);hold on % create 20 plots in 4 by 5 arrangement
%     for bb = 1:(length(my_Bin)-1) % loop over 1 to 11 because we want 12 interval from 13 numbers, 0 and 360 is the same here
% 
%         thisidx = all_data.subj_all==ss & all_data.use_trial==1 & location_Bin==bb;
%         thisidx2 = find(thisidx);% minimize the data size for plot by selecting a group of data with the above criterions
% 
% 
%         plot((my_Bin(bb+1)+my_Bin(bb))/2,mean(all_data.s_all.f_sacc_err(thisidx2)),'o','Color',bin_colors(bb,:),'MarkerFaceColor',bin_colors(bb,:),'MarkerSize',4);
%     end %plot(data point on x axis is ploted on the center of each interval e.g., (0+30)/2 = 15 and so on for all bin
% 
%     text(330,4.5,u_subj{ss},'FontAngle','italic','HorizontalAlignment','right');
% 
%     xlim([0 360]);ylim([0 5]);
%     xticks([0 90 180 270 360]);
%     yticks([0 10]);
%     % set(gca,'TickDir','out','XTick',[-180 0 180],'YTick',[-180 0 180]);
%     if ss == 1; xlabel('Target location (°)'); ylabel('Final Saccade Error (MAE; °)'); end
%     set(gca,'TickDir','out','LineWidth',0.75);
% 
% end
% 
% %same graph collapse across sub
% figure;
% for ss = 1:length(u_subj) 
%     for bb = 1:(length(my_Bin)-1) 
%         subplot(1,1,1);hold on
%         thisidx = all_data.subj_all==ss & all_data.use_trial==1 & location_Bin==bb;
%         thisidx2 = find(thisidx);% minimize the data size for plot by selecting a group of data with the above criterions
% 
%         plot((my_Bin(bb+1)+my_Bin(bb))/2,mean(all_data.s_all.f_sacc_err(thisidx2)),'o','Color',bin_colors(bb,:),'MarkerFaceColor',bin_colors(bb,:));
%         
%     end 
% 
%     xlim([0 360]);ylim([0 10]);
%     xticks([0 90 180 270 360]);
% 
%     xlabel('Target location (°)'); ylabel('Final Saccade Error (MAE; °)');
%     set(gca,'TickDir','out');
% end

%% plot targdiff against error

for ss= 1:length(u_subj)
    subplot(1,length(u_subj),ss); hold on
%     plot(targ_diff(all_data.subj_all == ss), all_data.s_all.f_sacc_align180(all_data.subj_all == ss),'ko','MarkerFaceColor','k','');
    scatter(targ_diff(all_data.subj_all == ss & all_data.use_trial==1), all_data.s_all.f_sacc_align180(all_data.subj_all == ss & all_data.use_trial==1),10,'k','filled','MarkerFaceAlpha',0.25);%MarkerFaceColor','k','');
    if ss ==1
    xlabel('target difference between current trial and previous trial in angle')
    ylabel('final saccade error in angle with algined target location')
    
    end
    
    title(u_subj{ss})
    ylim([-15,15]);
    axis square;
end
%% practice SD shading graph
my_Bin = -180:10:180;
diff_Bin = discretize(targ_diff,my_Bin);%put target location data into bin from 1-4, e.g., 1 is 0-90 degree and so on.

for ss= 1:length(u_subj)

    subplot(3,7,ss); hold on;

    for bb = 1:(length(my_Bin)-1)

        thisidx = all_data.subj_all==ss & all_data.use_trial==1 & diff_Bin==bb;

%         plot((my_Bin(bb+1)+my_Bin(bb))/2,mean(all_data.s_all.f_sacc_align180(thisidx)),'ko');%,'Color',bin_colors(bb,:),'MarkerFaceColor',bin_colors(bb,:),'MarkerSize',4);

        y = mean(all_data.s_all.f_sacc_align180(thisidx)); % your mean vector;
        x =(my_Bin(bb+1)+my_Bin(bb))/2;
        SEM = std(all_data.s_all.f_sacc_align180(thisidx))/sqrt((length(all_data.s_all.f_sacc_align180(thisidx))));
        curve1 = y + SEM;
        curve2 = y - SEM;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        fill(x2, inBetween, 'g');
        hold on;
        plot(x, y, 'r', 'LineWidth', 2);

        %SEM = std(all_data.s_all.f_sacc_align180(thisidx))/sqrt(length(all_data.s_all.f_sacc_align180(thisidx)).....
        if ss ==1
            xlabel('target difference between current trial and previous trial in angle')
            ylabel('final saccade error with algined target location')

        end
    end

    title(u_subj{ss})
    xlim([-180 180]);ylim([-10 10]);
    %     axis square;
    %set(gca,'TickDir','out','XTick',[-180 0 180],'YTick',[-180 0 180]);
end


%% SD plot with response bias against target location differences bin
%%% TODO: plot group moving average of response errors with standard error in shading  

my_Bin = -180:10:180; % create bin limit for -180 to 180, make 1 bin for 10 degree interval therefore result in 37 bins in total

%[BINS, EDGES] = discretize(X,N) where N is a scalar integer, divides the range of X into N uniform bins, and also returns the bin edges.
diff_Bin = discretize(targ_diff,my_Bin);%put target location data into bin from 1-4, e.g., 1 is 0-90 degree and so on.
bin_colors = hsv(length(my_Bin));%give each bin a distinct color

% which_subj_traj = 1; % to plot all, do which_subj_traj = u_subj;

for ss= 1:length(u_subj)

    subplot(3,7,ss); hold on;

    for bb = 1:(length(my_Bin)-1)

        thisidx = all_data.subj_all==ss & all_data.use_trial==1 & diff_Bin==bb;

        plot((my_Bin(bb+1)+my_Bin(bb))/2,mean(all_data.s_all.f_sacc_align180(thisidx)),'ko');%,'Color',bin_colors(bb,:),'MarkerFaceColor',bin_colors(bb,:),'MarkerSize',4);

%SEM = std(all_data.s_all.f_sacc_align180(thisidx))/sqrt(length(all_data.s_all.f_sacc_align180(thisidx)).....
        if ss ==1
            xlabel('target difference between current trial and previous trial in angle')
            ylabel('final saccade error with algined target location')

        end
    end

    title(u_subj{ss})
    xlim([-180 180]);ylim([-10 10]);
    %     axis square;
    %set(gca,'TickDir','out','XTick',[-180 0 180],'YTick',[-180 0 180]);
end
%% plot RT distribution and targdiff against RT

figure;
for ss = 1:length(u_subj)
    subplot(3,7,ss); hold on
    histogram(all_data.s_all.i_sacc_rt, 45,'BinLimits',[0 1.5]);

    xlabel('RT (s)');
    title(u_subj{ss});

end


figure;
for ss= 1:length(u_subj)
    subplot(1,length(u_subj),ss); hold on

    scatter(targ_diff(all_data.subj_all == ss & all_data.use_trial==1), all_data.s_all.i_sacc_rt(all_data.subj_all == ss & all_data.use_trial==1),10,'k','filled','MarkerFaceAlpha',0.25);%MarkerFaceColor','k','');
    if ss ==1
    xlabel('target difference between current trial and previous trial in angle')
    ylabel('RT (s)')
    
    end
    
    title(u_subj{ss})
%     ylim([-15,15]);
%     axis square;
end


%% plot individual sub data
% this_subj = 19;
% figure;
%  scatter(targ_diff(all_data.subj_all == this_subj & all_data.use_trial==1), all_data.s_all.f_sacc_align180(all_data.subj_all == this_subj & all_data.use_trial==1),10,'k','filled','MarkerFaceAlpha',0.25);%MarkerFaceColor','k','');
% 
% hold on
% 
% c = sqrt(2)/exp(-0.5);
% a = 3;
% w = 0.015;
% 
% x = [-180:1:180];
% 
% y = x.*a.*w.*c.*exp(-(w.*x).^2);
% 
% xlim([-180 180]);ylim([-10 10]);
% 
% plot(x,y,'r-');

%% manual error measurement

opt_param = nan(length(u_subj),2);
opt_ssRes = nan(length(u_subj),1);
opt_ve = nan(length(u_subj),1);

for ss = 1:length(u_subj)
 subplot(3,7,ss); hold on
% this_subj = 19;
my_x = targ_diff(all_data.subj_all == ss & all_data.use_trial==1);
my_y = all_data.s_all.f_sacc_align180(all_data.subj_all == ss & all_data.use_trial==1);

fun_x = my_x;
% fun_y = fun_x.*a.*w.*c.*exp(-(w.*fun_x).^2);
plot (my_x,my_y,'r.');hold on
% % plot (fun_x,fun_y,'b.');
% 
% error_mea= my_y - fun_y;
% square_error = error_mea.^2;
% sum_error = sum(square_error);

% w is p(1); a is p(2)
myerrfcn = @(p) sum( (  my_y -  (fun_x.*p(2).*p(1).*c.*exp(-(p(1).*fun_x).^2)) ).^2  ) ; % computes SSE

options = optimset('fmincon');
options.Display = 'off';

% init_guess = [1/15 10];  
init_guess = [0.015 3]; %0.015 3;

[bf_params, bf_err, ex_flag] = fmincon( myerrfcn, init_guess, [], [], [], [], [0 -inf], [inf inf], [], options );
%[bf_params, bf_err, ex_flag] = fmincon( myerrfcn, init_guess);

opt_param(ss,:) = [bf_params(1), bf_params(2)]; 
opt_ssRes(ss) = [bf_err];
opt_ve(ss) = [1-opt_ssRes(ss)/sum((my_y-mean(my_y)).^2)];

c = sqrt(2)/exp(-0.5);
a = bf_params(2);
w = bf_params(1);

x = [-180:1:180];

y = x.*a.*w.*c.*exp(-(w.*x).^2);

xlim([-180 180]);ylim([-10 10]);
plot(x,y,'b-');
    if ss ==1
    xlabel('target difference (°)')
    ylabel('final saccade error (°)')
    
    end
    
    title(u_subj{ss})
    ylim([-10,10]);
    axis square;
clear a w y bf_params;

end


%% parameter fit

good_sub = opt_ve >= 0;% excluded two subj (16 and 17 (sub039,sub043) because the model fit explain the variance less than the data itself

figure;

subplot(1,2,1);hold on

plot(0, opt_param(good_sub,1),'ko');
subplot(1,2,2);
plot(1, opt_param(good_sub,2),'ro')
figure;
scatter(opt_param(good_sub,1),opt_param(good_sub,2))

%% loop over to find better initial guess
% guess range for a [0 0.032 8], w [0 0.004 1];
% genenrate two sets of 250 numbers in the above two ranges, make random combination of the two which would be 62500 set of combination.
% then loop over them in initial guess, and generate an opt_ve value every
% loop to find the biggest average opt_ve across subj
% therefore, three numbers needed every loop: the initial guess and the
% average opt_ve


%% plot SD amplitude by subj

figure;
for ss = 1:length(u_subj)
    
    plot(good_subNum.',opt_param(:,2),'ko')
    xlabel('subjects')
    ylabel('Serial Dependence (°)')
  
    ylim([-2,8]);
end
% %% plot error against ITI
% 
% iti_all = all_data.c_all(:,2);
% figure;
% for ss= 1:length(u_subj)
%     subplot(3,7,ss); hold on
% 
%     scatter(iti_all(all_data.subj_all == ss & all_data.use_trial==1), all_data.s_all.f_sacc_align180(all_data.subj_all == ss & all_data.use_trial==1),10,'k','filled','MarkerFaceAlpha',0.25);%MarkerFaceColor','k','');
% 
%     xlabel('ITI')
%     ylabel('Final Saccade Error')
% 
% 
%     title(u_subj{ss})
%      xlim([1 5]);ylim([0  10]);
%     %     axis square;
% end
