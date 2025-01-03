% Load preprocessed EEG set and extract FRP and features without overlap
% correction
% RVS 2024

clear all; close all force
init_unfold

% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
%%%%%%%%%
sublist = 19:181; % TODO: replace with full sublist, short list used for dev
%%%%%%%%%
exclude_linenoise = [30 36 98 101 102 109 111 114 118 122 125 131 134 136 139]; % TODO: deal with this line noise
exclude_noise = [32 86]; %n other noise such as excessive jumps, blinks -  131? , movement
exclude_movement = [];
exclude_missingevents = [52 57 73 111 120 153]; % 57 73 because no eyetracking, others because EEG stopped recording early
exclude_noEEG = [1:18 23 77 88 138 79 87 92 127 33 129 152]; % 79 onwards missing from reparsed fixations
exclude_other = [22:24 26 27 31 39 40 78 159 160 173 174 179 180]; % TODO find reason for these - no .set why? # 179 no gaze, 181 no reparsed fixations
exclude_glmfitfail = [19 35 68 147 149 ] ; % 54 onwards are only failing with sac splines
exclude = unique([exclude_linenoise exclude_noise exclude_movement exclude_missingevents exclude_noEEG exclude_glmfitfail exclude_other]); % Subj to exclude because no eeg or no trigger etc.
sublist = sublist(~ismember(sublist,exclude) );
%sublist = sublist( ismember(sublist,find(hasTriggerList.sdcard==1)));

%dir_raw = '/Volumes/Blue1TB/localEyeMindLink/Data';
dir_events = '~/Emotive Computing Dropbox/Rosy Southwell/EyeMindLink/Processed/events/';
dir_pre = fullfile('/Volumes/Blue1TB/EEG_processed') ; % prepro in MNE
dir_in = fullfile(dir_pre, 'preprocessed_set');
dir_fif = fullfile(dir_pre, 'preprocessed_fif');

unfdir = '/Volumes/Blue1TB/EEG_processed/unfolded_FRP_reparsed_v7/';
featdir = '/Volumes/Blue1TB/EEG_processed/unfolded_FRP_reparsed_v7/n400_stats_nodc/';
mkdir(featdir)

times = -300:10:790; % entire FRP time basis
cols = {'type','latency','duration','task','identifier','page_fixation_ix','IA_ID'};
FRPavg=[]; FRPall=[];
for s = 1:length(sublist)
        tic;
        clear EEG stats

        pID = ['EML1_',sprintf('%03d',sublist(s))];
        try
            % read unfstats because this has urevent for selecting trials
            % that passed QC
            stats = readtable(fullfile(unfdir, [pID '_reading_N400_stats.csv']));
        catch
            disp(['Failed to load stats for ' pID])
            continue
        end
        % events = dc_EEG.event;
        try
            % read EEG set
             EEG = pop_loadset(fullfile(unfdir, [pID '.set']));
        catch
            disp(['Failed to load EEG for ' pID])
            continue
        end

        % Epoch EEG around fixations
        EEG_epoch = uf_epoch(EEG,'timelimits',[-0.3 0.8]);
        trials = stats.urevent;
        % select only events from urevent
        keep_evt = struct2table(EEG.event(trials,:));
        keep_evt.urevent=trials;

        nodc_FRP = EEG_epoch.data(:,:,trials);
        stats_old = stats;

        stats_nodc = stats(:,cols);
        
        %% loop over channels and recompute features per fixation
        channels = 1:8;
        chanlabels = {'CPz', 'FCz'  , 'AFF5h' ,'AFF6h' , 'CCP5h' ,'CCP6h' , 'PPO9h' , 'PPO10h'};

        for c = channels
            chanlabel = chanlabels{c};

            % N400 features
            win=[300 500];
            sel = find(times<= win(2) & times>=win(1));
            y=squeeze(nodc_FRP(c,sel,:));

            % absolute max magnitude (regardless of polarity)
            [mag, ix] = max(abs(y),[],1);
            mag = zeros(size(ix(:))); % max val for each trial. I'm being dumb but cant work out how to vectorise this rn
            lat = squeeze(times(sel(squeeze(ix))));
            keep_evt.(['n400_magnitude_' chanlabel]) = mag;
            keep_evt.(['n400_latency_' chanlabel]) = lat';

            % max and min magnitude
            [mag, ix] = max(y,[],1);
            lat = squeeze(times(sel(squeeze(ix))));
            stats.(['n400_max_magnitude_' chanlabel]) = mag';
            stats.(['n400_max_latency_' chanlabel]) = lat';

            [m, ix] = min(y,[],1);
            mag = zeros(size(ix(:))); % min val for each trial. I'm being dumb but cant work out how to vectorise this rn
            for z = 1:length(ix(:))% loop over trials
                mag(z) = y(ix(z),z);
            end
            lat = squeeze(times(sel(squeeze(ix))));
            stats.(['n400_min_magnitude_' chanlabel]) = mag;
            stats.(['n400_min_latency_' chanlabel]) = lat';
            
            % zero crossings
            [rate, zc] = zerocrossrate(y,Method="comparison");
            stats.(['n400_zero_crossings_' chanlabel]) = zc';

            % mean value
            stats.(['n400_mean_' chanlabel]) = mean(y,1)';

            %% P1 features (for latency see https://www.mdpi.com/2411-5150/4/1/11)
            win=[70 120];
            sel = find(times<= win(2) & times>=win(1));
            y=squeeze(nodc_FRP(c,sel,:));
            % mean value
            stats.(['p1_mean_' chanlabel]) = mean(y,1)';

            %% N1 features (for latency see https://www.mdpi.com/2411-5150/4/1/11)
            win=[140 280]; 
            sel = find(times<= win(2) & times>=win(1));
            y=squeeze(nodc_FRP(c,sel,:));
            % mean value
            stats.(['n1_mean_' chanlabel]) = mean(y,1)';

        end

        %% get average FRP 
        FRPavg=mean(nodc_FRP,3);
        FRPall(:,:,s) =FRPavg;
    
        writetable(stats, fullfile(featdir,[pID '_reading_N400_stats.csv'])); 

end
colororder("glow12")
figure(1)
plot(times,mean(FRPall,3),'LineWidth',2);
hold on
plot(times, zeros(size(times)))
legend(chanlabels)

