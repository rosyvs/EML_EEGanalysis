% Load EEG and align it to eyetracker data
% perform eyetracker-optimised ICA (OPTICAT)
%------------------------------------------

% eyetracker timestamps are milliseconds measured from the time when the tracker software was started
% eeg timestamps are sample number but the sampling rate is 1000Hz here so
% no conversion necessary. **warning** IF reusing this script, you will need to be
% cautious about wehther you are using time in sample or in milliseconds
% and convert where appropriate.

% Load preprocessed data from MNE-python
% - this data has been cleaned using robust detrending, DSS_based line
%   noise removal, and has events labelled from the trial log file.
% - ** It might be worth porting all the MNE-python prepro back to Matlab
%   seeing as it mostly relies on libraries adapted from Matlab in the first
%   place **

clear all; close all
eeglab nogui % sets path defaults
addpath(genpath('C:\Users\roso8920\Documents\MATLAB\eeglab_current\eeglab2021.0\plugins\eye-eeg0.85'))
% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
sublist = find(hasTriggerList.sdcard==1);
sublist = sublist(sublist~=73); % no eyetracking for these subjects
eeg_exclude = [20, 21,22, 26,77]; % Subj to exclude because no eeg or no trigger

dir_raw = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EyeMindLink\Data\';
dir_pre = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EML Rosy\Data\EEG_processed';
mkdir( dir_pre, 'opticat_cleaned')
for s = 1:length(sublist)
    close all
    tic;
    pID = ['EML1_',sprintf('%03d',sublist(s))];    
    
    % use Fieldtrip to load fif data from MNE (TODO can EEGlab read fif
    % directly??)
    cfg = [];
    cfg.dataset = fullfile(dir_pre, [pID '_p.fif']);
    EEG_ft = ft_preprocessing(cfg);
    
    % Read info txt to determine whether EEG+triggers are from SD card
    % (default) or streamed (backup)
    triginfo = readtxtfile(fullfile(dir_pre, [pID '-info.txt']));
    
    % read events.csv for triggers and descriptions
    logtrig = readtable(fullfile(dir_raw,pID,'EEG','events.csv'));
    % copy the correct EEGsample column for use depending on triginfo
    if contains(triginfo, 'LA0','IgnoreCase',false)
        logtrig.eeg_use_sample = logtrig.eegSD_sample_est;
    else
        logtrig.eeg_use_sample = logtrig.eeg_sample_est;
    end
    %% read eyetracker events
    eye_srate=1000; % Hz sampling rate of eyetracker data
    blinkfiles = dir(fullfile(dir_raw,pID,'Unpacked','*Blink.csv'));
    blinks=[];
    for i = 1:length(blinkfiles)
        blinks =  [blinks;...
            readtable(fullfile(dir_raw,pID,'Unpacked',blinkfiles(i).name))];
    end
    
    saccfiles = dir(fullfile(dir_raw,pID,'Unpacked','*Saccade.csv'));
    saccades=[];
    for i = 1:length(saccfiles)
        saccades =  [saccades;...
            readtable(fullfile(dir_raw,pID,'Unpacked',saccfiles(i).name))];
    end
    
    fixfiles = dir(fullfile(dir_raw,pID,'Unpacked','*Fixation.csv'));
    fixations = [];
    for i = 1:length(fixfiles)
        fixations = [fixations;...
            readtable(fullfile(dir_raw,pID,'Unpacked',fixfiles(i).name))];
    end
    % drop extreme fixations
    fixations = fixations(fixations.duration>100 & fixations.duration<1000,:);
    
    
    msgfiles = dir(fullfile(dir_raw,pID,'Unpacked','*Message.csv'));
    messages=[];
    for i = 1:length(msgfiles)
        messages =  [messages;...
            readtable(fullfile(dir_raw,pID,'Unpacked',msgfiles(i).name))];
    end
    
    %% get eyetracker timestamp rel to eeg timestamp by matching log to ET messages
    % warning: this is messed up for subjects  33, 55, 56, 57, 58, 70
    % because the coutner went back to 0 for eyelink samples. ask Megan for
    % fixed. 
    
    messages_trialonsets = messages(contains(messages.text,'TRIALID ') & ~contains(messages.text,{'DriftCorrect','Recal'}),:);
    %logtrig_trialonsets = logtrig(~contains(logtrig.EVENT, {'DriftCorrect','Recal'}),:)
    logtrig.eye_sample = NaN(height(logtrig),1);
    % for each eyetracker message, find its corresponding row in logtrig
    for i = 1:height(messages_trialonsets)
        event = messages_trialonsets.text{i}; event=event(9:end); % remove 'TRIALID '
        log_ix = find(matches(logtrig.EVENT,event));
        eye_sample = messages_trialonsets.time(i) ;
        temp(i) = length(log_ix);

        if temp(i) > 1
            warning(['multiple log entries found for eyetracker message ' event]);
            % TODO find a way to sensibly choose in case of duplicates
           %  log_ix = log_ix(logtrig.eye_sample() > messages_trialonsets(i-1)); 
           log_ix=log_ix(end); % most likely the exp was started falsely once and the earlier one(s) should be discarded
        end
        if isempty(log_ix)
            warning(['No log event matches eyetracker message: ' event])
        else
            logtrig.eye_sample(log_ix) = eye_sample;

        end
      
    end
    % check eyetracker vs log timing jitter
    logtrig.eye_diff_since_last = [0; diff(logtrig.eye_sample)];
    eye_log_jitter = range(logtrig.diff_since_last - logtrig.eye_diff_since_last);
    % check EEG vs log timing jitter - remember we take the SD card
    % recording over the streamed recording
    logtrig.eeg_use_diff_since_last = [0; diff(logtrig.eeg_use_sample)];
    eeg_log_jitter = range(logtrig.diff_since_last - logtrig.eeg_use_diff_since_last);

    
    
    %% task events into event structure
    % type, latency,duration
    task_events = logtrig(:,{'EVENT','eeg_use_sample'});
    task_events = renamevars(task_events,{'EVENT','eeg_use_sample'},{'type','latency'});
   
    %% interpolate all eye event timestamps to get closest eeg sample
    % first select rows with non-NaN values for both eeg and eyetracker samples
    ix = ~isnan(logtrig.eye_sample + logtrig.eeg_use_sample);
    F = griddedInterpolant(logtrig.eye_sample(ix), logtrig.eeg_use_sample(ix), 'spline');
    % plot interpolant: this should be a straight line and the red points shold
    % be on the line
    h0=figure()
    plot(min(F.GridVectors{1,1}):max(F.GridVectors{1,1}),F(min(F.GridVectors{1,1}):max(F.GridVectors{1,1})))
    hold on
    plot(logtrig.eye_sample(ix), logtrig.eeg_use_sample(ix),'ro')
    xlabel('eyetracker sample')
    ylabel('EEG sample')
    
    saveas(h0, fullfile(dir_pre, 'opticat_cleaned', [pID, '_eeg-eye_sample_interpolant.png']))
    
    blinks.eeg_tStart = round(F(blinks.tStart));
    saccades.eeg_tStart = round(F(saccades.tStart));
    fixations.eeg_tStart = round(F(fixations.tStart));
    
     %%%% TODO
    % choice of eye to use - earliest event? both? best tracked eye?
    
    
    %% convert to EEGLAB to use OPTICAT
    EEG = fieldtrip2eeglab(EEG_ft.hdr, EEG_ft.trial{1});
    EEG = pop_select(EEG,'nochannel',{'x_dir','y_dir','z_dir'});
    % note data units appears to be V here. EEGlab expects microvolt
    EEG.data = EEG.data*10e6; % convert V to uV
    
    EEG_copy = EEG; % make a copy of the pre-cleaned data for comparison
    
    %% OPTICAT ICA
    % This script implements the procedures from:
    % Dimigen, O. (2020). Optimizing the ICA-based removal of ocular artifacts
    % from free viewing EEG. NeuroImage, https://doi.org/10.1016/j.neuroimage.2019.116117
    
    % PARAMETERS:
    opticat_hp = 2; % highpass filter in Hz used for ICA training
    opticat_lp = 100;
    opticat_downsample = 200; % downsample to rate (Hz)
    opticat_epoch_dur = 3; % arbitrary epochs to cut the data into (because we don't have a straightforward event-related dataset)
    opticat_reject_uV = 500; % reject pseudoepochs from training data with range exceeding this threshold
    opticat_overweighting = 1; % additional spike potentials as proportion of orig data
    opticat_mean_centre = 'true';
    eeg_channels = 1:8; % use all 8 EEG channels - the other 3 are accelerometer
    
    
    %% Add eye events to EEGlab event structure
    blinks.type =  repmat({'blink_either_eye'},height(blinks),1);
    blinks.latency = blinks.eeg_tStart;
    
    fixations.type = repmat({'fixation_either_eye'},height(fixations),1);
    fixations.latency = fixations.eeg_tStart;
    
    saccades.type = repmat({'saccade_either_eye'},height(saccades),1);
    saccades.latency = saccades.eeg_tStart;
    
    ev_fields = {'type','latency','duration'};
    eye_events = [blinks(:,ev_fields); fixations(:,ev_fields); saccades(:,ev_fields)];
    eye_events = sortrows(eye_events,'latency');
    % truncate at end of eeg recording
    eye_events = eye_events(eye_events.latency+eye_events.duration <= EEG.pnts,:);
    
    EEG.event = table2struct(eye_events);
    EEG_copy.event = table2struct(eye_events);
    
    %% Make pseudoepoch event structure
    n_pseudo = floor(EEG.pnts/(opticat_epoch_dur*EEG.srate));
    
    event_pseudo = table();
    event_pseudo.latency =  [0:n_pseudo-1]'*opticat_epoch_dur*EEG.srate+1;
    event_pseudo.type =repmat({'pseudo_epoch'}, height(event_pseudo),1);
    event_pseudo.duration = repmat(opticat_epoch_dur*EEG.srate, height(event_pseudo),1);
    
    event_training = [eye_events; event_pseudo];
    event_training = sortrows(event_training, {'latency'});
    
    EEG_training = EEG;
    EEG_training.event = table2struct(event_training);
    
    
    
    %% Prepare training data
    % filter data copy for training
    EEG_training = pop_eegfiltnew(EEG_training,opticat_hp,opticat_lp);
    % resample data copy for training
    EEG_training = pop_resample(EEG_training, opticat_downsample);
    % check event structure
    EEG_training = eeg_checkset(EEG_training);
    % apply pseudoepoching
    EEG_training = pop_epoch(EEG_training, {'pseudo_epoch'}, [0 opticat_epoch_dur]);
    %% Reject +- 500 uV pseudoepochs
    EEG_training = pop_eegthresh(EEG_training, 1, eeg_channels, -500, 500, 0, opticat_epoch_dur, 1,0);
    % eye events in the event field are already numbered with corresponding
    % pseudoepoch number, which means corresponding eye
    % data will also be removed from training set
    % TODO count trials before and after
    reject_events = ismember([EEG_training.event(:).epoch], find(EEG_training.reject.rejthresh)) ;
    EEG_training = pop_rejepoch(EEG_training,find(EEG_training.reject.rejthresh) ,0);
    
    
    
    %% Overweight spike potentials
    %   Repeatedly append intervals around saccade onsets (-20 to +10 ms) to training data
    EEG_training = pop_overweightevents(EEG_training,'saccade_either_eye',[-0.02 0.01],opticat_overweighting,opticat_mean_centre);
    
    %% Run ICA
    fprintf('\nRunning ICA on optimized training data...')
    EEG_training = pop_runica(EEG_training,'extended',1,'interupt','on','chanind',eeg_channels); % or use binary ICA for more speed
    
    % Remember ICA weights & sphering matrix
    wts = EEG_training.icaweights;
    sph = EEG_training.icasphere;
    
    %% Transfer unmixing weights
    % Remove any existing ICA solutions from your original dataset
    EEG.icaact      = [];
    EEG.icasphere   = [];
    EEG.icaweights  = [];
    EEG.icachansind = [];
    EEG.icawinv     = [];
    
    fprintf('\nTransfering ICA weights from training data to original data...')
    EEG.icasphere   = sph;
    EEG.icaweights  = wts;
    EEG.icachansind = eeg_channels;
    EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv
    
    %% Eye-tracker-guided selection of ICs
    fprintf('\nIdentifying ocular ICs via saccade/fixation variance-ratio threshold...')
    
    IC_THRESHOLD     = 1.1;   % variance ratio threshold (determined as suitable in Dimigen, 2020)
    SACC_WINDOW      = [5 0]; % saccade window (in samples!) to compute variance ratios (see Dimigen, 2020)
    PLOTFIG          = true;  % plot a figure visualizing influence of threshold setting?
    ICPLOTMODE       = 0;     % plot component topographies (inverse weights)? (2 = only plot "bad" ocular ICs)
    FLAGMODE         = 3;     % overwrite existing rejection flags? (3 = yes)
    MAXCOMP = 2; % how many components to remove MAXIMUM? 
    
    % Automatically flag ocular ICs (Plöchl et al., 2012)
    [EEG, varratiotable] = pop_eyetrackerica(EEG,'saccade_either_eye','fixation_either_eye',SACC_WINDOW,IC_THRESHOLD,FLAGMODE,PLOTFIG,ICPLOTMODE);
    h=gcf();
    saveas(h, fullfile(dir_pre, 'opticat_cleaned', [pID, '_varratios.png']))

    % Remove flagged ocular ICs
    badcomps = EEG.reject.gcompreject;
    % choose max 2 with highest variance
    [vr,worst2comps] = maxk(varratiotable(:,3),MAXCOMP);
    badcomps = intersect(find(badcomps), worst2comps);
    
    % plot chosen components & remaining components    
    h_compfig=figure(22); clf
    for ci=1:length(badcomps)
    subplot(1,MAXCOMP,ci)
    topoplot(EEG.icawinv(:,badcomps(ci)),EEG.chanlocs, 'verbose', ...
			      'off', 'chaninfo', EEG.chaninfo, 'numcontour', 8);
    title(['Component ' num2str(badcomps(ci)) '| Var ratio = ' num2str(varratiotable(badcomps(ci),3))])
    end
    sgtitle([ pID ' eyeca removed components'],'Interpreter','none')
    saveas(h_compfig, fullfile(dir_pre, 'opticat_cleaned', [pID, '_detected_components.png']))

    %% remove bad components
    EEG      = pop_subcomp(EEG,badcomps); % remove them    
    
    %%  Plot SRP + FRP before and after artefact removal
    % Compute saccade-locked ERPs from clean data, then baseline (ugh mixed units - seconds then milliseconds)
    epoch_sac_clean = pop_rmbase(pop_epoch(EEG,{'saccade_either_eye'},[-0.2 0.6]),[-200 0]);
    epoch_sac = pop_rmbase(pop_epoch(EEG_copy,{'saccade_either_eye'},[-0.2 0.6]),[-200 0]);
    epoch_fix_clean = pop_rmbase(pop_epoch(EEG,{'fixation_either_eye'},[-0.2 0.6]),[-200 0]);
    epoch_fix = pop_rmbase(pop_epoch(EEG_copy,{'fixation_either_eye'},[-0.2 0.6]),[-200 0]);
    
    h1=figure('name','Saccade-related potentials before and after OPTICAT-ICA');
    subplot(1,2,1)
    plot(epoch_sac_clean.times, mean(epoch_sac_clean.data(eeg_channels,:,:),3))
    ylabel('Saccade-related ERP')
    xlabel('Time after saccade [ms]')
    title([pID ': orig data' ], 'Interpreter','none')
    subplot(1,2,2)
    plot(epoch_sac.times, mean(epoch_sac.data(eeg_channels,:,:),3))
    ylabel('Saccade-related ERP')
    xlabel('Time after saccade [ms]')
    title(['clean data - ' num2str(length(badcomps)) ' removed components'])
    legend({EEG.chanlocs.labels})
    set(h1,'Position',[190 340 1200 400])
    
    h2=figure('name','Fixation-related potentials before and after OPTICAT-ICA');
    subplot(1,2,1)
    plot(epoch_fix_clean.times, mean(epoch_fix_clean.data(eeg_channels,:,:),3))
    ylabel('Fixation-related ERP')
    xlabel('Time after fixation [ms]')
    title([pID ': orig data' ], 'Interpreter','none')
    subplot(1,2,2)
    plot(epoch_fix.times, mean(epoch_fix.data(eeg_channels,:,:),3))
    ylabel('Fixation-related ERP')
    xlabel('Time after fixation [ms]')
    title(['clean data - ' num2str(length(badcomps)) ' removed components'])
    legend({EEG.chanlocs.labels})
    set(h2,'Position',[190 340 1200 400])
    
    %% save data and plots
    saveas(h1, fullfile(dir_pre, 'opticat_cleaned', [pID, '_saccades.png']))
    saveas(h2, fullfile(dir_pre, 'opticat_cleaned', [pID, '_fixations.png']))
    pop_saveset(EEG, 'filename',pID, 'filepath',fullfile(dir_pre, 'opticat_cleaned'), 'savemode','onefile')
    disp(['Done EyeCA for ' pID])
    
    % Finished ICa for 1 subject! 
    load train 
    sound(y,Fs)

    toc
end