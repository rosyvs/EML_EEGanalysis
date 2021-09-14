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
% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
sublist = [60:100];
exclude = [20:26, 54, 73, 77]; % Subj to exclude because no eeg or no trigger etc. 54 i fixable - check message file 6
sublist = sublist(~ismember(sublist,exclude) & ismember(sublist,find(hasTriggerList.sdcard==1)));

dir_raw = '/Volumes/Blue1TB/EyeMindLink/Data';
dir_pre = fullfile('..','..','Data','EEG_processed') ;
dir_fif = fullfile(dir_pre, 'preprocessed_fif');
mkdir( dir_pre, 'preprocessed_set')

for s = 1:length(sublist)
    close all
    tic;
    pID = ['EML1_',sprintf('%03d',sublist(s))];    
    
    % use Fieldtrip to load fif data from MNE (TODO can EEGlab read fif
    % directly??)
    cfg = [];
    cfg.dataset = fullfile(dir_fif, [pID '_p.fif']);
    EEG_ft = ft_preprocessing(cfg);
        %% convert to EEGLAB 
    EEG = fieldtrip2eeglab(EEG_ft.hdr, EEG_ft.trial{1});
    EEG = pop_select(EEG,'nochannel',{'x_dir','y_dir','z_dir'});
    % note data units appears to be V here. EEGlab expects microvolt
    EEG.data = EEG.data*10e6; % convert V to uV
    
    %% info
    % Read info txt to determine whether EEG+triggers are from SD card
    % (default) or streamed (backup)
    triginfo = readtxtfile(fullfile(dir_fif, [pID '-info.txt']));
    
    % read events.csv for triggers and descriptions
    logtrig = readtable(fullfile(dir_raw,pID,[pID '_events.csv']));
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
  %  fixations = fixations(fixations.duration>100 & fixations.duration<1000,:);
    
    
    msgfiles = dir(fullfile(dir_raw,pID,'Unpacked','*Message.csv'));
    messages=[];
    for i = 1:length(msgfiles)
        messages =  [messages;...
            readtable(fullfile(dir_raw,pID,'Unpacked',msgfiles(i).name))];
    end
        
    %% task events into event structure
    % type, latency,duration
    task_events = logtrig(:,{'EVENT_logtrig','eeg_use_sample','duration_sec'});
    task_events = renamevars(task_events,{'EVENT_logtrig','eeg_use_sample','duration_sec'},{'type','latency','duration'});
   task_events.duration = task_events.duration * EEG_ft.fsample; % into samples
   
    %% interpolate all eye event timestamps to get closest eeg sample
    % first select rows with non-NaN values for both eeg and eyetracker samples
    ix = ~isnan(logtrig.eye_sample + logtrig.eeg_use_sample);
    F = griddedInterpolant(logtrig.eye_sample(ix), logtrig.eeg_use_sample(ix), 'spline');
     
    blinks.eeg_tStart = round(F(blinks.tStart));
    saccades.eeg_tStart = round(F(saccades.tStart));
    fixations.eeg_tStart = round(F(fixations.tStart));
    
     %%%% TODO
    % choice of eye to use - earliest event? both? best tracked eye?
    

    
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
    
    events = [eye_events; task_events];
    events=sortrows(events,'latency');
    EEG.event = table2struct(events);
    
    %% save .set
     pop_saveset(EEG, 'filename',pID, 'filepath',fullfile(dir_pre, 'preprocessed_set'), 'savemode','onefile')

    disp(['Converted .fif to .set for ' pID])
    
  
    toc
end