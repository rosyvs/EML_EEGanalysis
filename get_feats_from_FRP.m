% Load existing FRPs extracted from EEG using unfold
% extract new features
% RVS 2024

clear all; close all force

% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
%%%%%%%%%
sublist = 19:158; % TODO: replace with full sublist, short list used for dev
%%%%%%%%%
exclude_linenoise = [30 36 98 101 102 109 111 114 118 122 125 131 134 136 139]; % TODO: deal with this line noise
exclude_noise = [32 86]; %n other noise such as excessive jumps, blinks -  131? , movement
exclude_movement = [];
exclude_missingevents = [52 57 73 111 120 153]; % 57 73 because no eyetracking, others because EEG stopped recording early
exclude_noEEG = [1:18 23 77 88 138];
exclude_other = [22:24 26 27 31 39 40 78 160]; % TODO find reason for these - no .set why?
exclude_glmfitfail = [35 68 147 149];
exclude = unique([exclude_linenoise exclude_noise exclude_movement exclude_missingevents exclude_noEEG exclude_glmfitfail exclude_other]); % Subj to exclude because no eeg or no trigger etc.
sublist = sublist(~ismember(sublist,exclude) );
%sublist = sublist( ismember(sublist,find(hasTriggerList.sdcard==1)));

%dir_raw = '/Volumes/Blue1TB/localEyeMindLink/Data';
dir_events = '~/Emotive Computing Dropbox/Rosy Southwell/EyeMindLink/Processed/events/';
dir_pre = fullfile('/Volumes/Blue1TB/EEG_processed') ; % prepro in MNE
dir_in = fullfile(dir_pre, 'preprocessed_set');
dir_fif = fullfile(dir_pre, 'preprocessed_fif');

unfdir = '/Volumes/Blue1TB/EEG_processed/unfolded_FRP_reparsed_v5';
featdir = '/Volumes/Blue1TB/EEG_processed/unfolded_FRP_reparsed_v5/n400_stats_recomputed/';
mkdir(featdir)

times = -300:10:790; % entire FRP time basis

for s = 1:length(sublist)
        tic;

        pID = ['EML1_',sprintf('%03d',sublist(s))];
        try
        stats = readtable(fullfile(unfdir, [pID '_reading_N400_stats.csv']));
        catch
            disp(['Failed to load stats for ' pID])
            continue
        end
        % events = dc_EEG.event;
        

        channels = 1:8;
        chanlabels = {'CPz', 'FCz'  , 'AFF5h' ,'AFF6h' , 'CCP5h' ,'CCP6h' , 'PPO9h' , 'PPO10h'};
        % loop over channels and recompute features per fixation
        for c = channels
            chanlabel = chanlabels{c};
            dcFRP = readmatrix(fullfile(unfdir, [pID '_unfoldedFRP-reading_' chanlabel '.csv']));

     
            %% N400 features
            win=[300 500];
            sel = find(times<= win(2) & times>=win(1));
            y=dcFRP;
            y = y(sel,:);
            [m, ix] = max(y,[],1);
            mag = zeros(size(ix(:))); % max val for each trial. I'm being dumb but cant work out how to vectorise this rn
            for z = 1:length(ix(:))% loop over trials
                mag(z) = y(ix(z),z);
            end
            lat = squeeze(times(sel(squeeze(ix))));
            stats.(['n400_max_magnitude_' chanlabel]) = mag;
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

        end

        writetable(stats, fullfile(featdir,[pID '_reading_N400_stats.csv']));        %% Plot single trials
        
     
    

end

