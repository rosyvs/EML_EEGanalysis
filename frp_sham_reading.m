% fixations sham vs word

% cleaned vs uncleaned data

clear all; close all
eeglab nogui % sets path defaults
addpath(genpath('C:\Users\roso8920\Documents\MATLAB\eeglab_current\eeglab2021.0\plugins\eye-eeg0.85'))

% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
sublist = find(hasTriggerList.sdcard==1);
sublist = sublist(sublist~=73); % EEG but no eyetracking for these subjects

dir_raw = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EyeMindLink\Data\';
dir_pre = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EML Rosy\Data\EEG_processed\opticat_cleaned\';
dir_info = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EML Rosy\Data\EEG_processed\';

mkdir( dir_pre, 'FRPs_eyeCA')

%%
for s = 1:length(sublist)
    tic;
    close all
    pID = ['EML1_',sprintf('%03d',sublist(s))];
    EEGica = pop_biosig(fullfile(dir_pre, [pID '.edf']));
    
    
    
    %% TODO rerun ICA with save dataset as EEGLab %%
    
    
    
    %% read log
    % Read info txt to determine whether EEG+triggers are from SD card
    % (default) or streamed (backup)
    triginfo = readtxtfile(fullfile(dir_info, [pID '-info.txt']));
    
    % read events.csv for triggers and descriptions
    logtrig = readtable(fullfile(dir_raw,pID,'EEG','events.csv'));
    % copy the correct EEGsample column for use depending on triginfo
    if contains(triginfo, 'LA0','IgnoreCase',false)
        logtrig.eeg_use_sample = logtrig.eegSD_sample_est;
    else
        logtrig.eeg_use_sample = logtrig.eeg_sample_est;
    end
    
    
    %% fix events field
    % event types have been overwritten with uninformative names
    % 1 = blink, 2=fixation, 3=saccade
    oldstr = {'condition 1','condition 2','condition 3'};
    replacements = {'blink_either_eye','fixation_either_eye','saccade_either_eye'};
    events = struct2table(EEGica.event);
    events.type = replace(events.type,oldstr,replacements);
    events = removevars(events,'edftype');
    
    %% make new event structure with fixations by reading vs sham
    task_events = logtrig(:,{'EVENT','VAL','eeg_use_sample'});
    task_events = renamevars(task_events,{'EVENT','eeg_use_sample'},{'type','latency'});
    task_events.duration = [diff(task_events.latency);NaN];
    task_events = task_events(task_events.VAL==7|task_events.VAL==20,:);
    
    task_events.type(task_events.VAL==7) = {'reading'};
    task_events.type(task_events.VAL==20) = {'sham'};
    task_events.urevent=NaN(height(task_events),1);
    task_events = removevars(task_events,'VAL');
    
    for e = 1:height(task_events)
        sel = find((events.latency >= task_events.latency(e)) & (events.latency < task_events.latency(e)+task_events.duration(e) ));
        events.type(sel) = strcat(events.type(sel), '_WITHIN_', repmat(task_events.type(e),length(sel),1));
    end
    
    events = sortrows([events; task_events(:,{'type','latency','duration','urevent'})],'latency');
    
    % remove NaN latency events
    events = events(~isnan(events.latency),:)
    
    EEGica.event = table2struct(events);
    EEGica = eeg_checkset(EEGica);
    
    %% epoch into sham/reading trials
    % Preprocess data
    EEGp = pop_eegfiltnew(EEGica,0.1,100);
    EEGp = pop_resample(EEGp, 200);
    EEGp = eeg_checkset(EEGp);
    EEG_epoched = pop_epoch(EEGp, {'fixation_either_eye_WITHIN_reading','fixation_either_eye_WITHIN_sham'},[-0.2,0.8]);
    
    % subsample to equate n epochs
    
    %% Plot epochs
    % pop_erpimage(EEG_epoched,1)
    
    epoch_reading = pop_rmbase(pop_epoch(EEGp,{'fixation_either_eye_WITHIN_reading'},[-0.2 0.8]),[-100 0]);
    epoch_sham = pop_rmbase(pop_epoch(EEGp,{'fixation_either_eye_WITHIN_sham'},[-0.2 0.8]),[-100 0]);
    
    
    h2=figure('name','Fixation-related potentials by event type');clf
    ax(1)=subplot(2,1,1);
    plot(epoch_reading.times, mean(epoch_reading.data,3))
    ylabel('Fixation-related ERP')
    xlabel('Time after fixation [ms]')
    title([pID ': Reading' ], 'Interpreter','none')
    ax(2)=subplot(2,1,2);
    plot(epoch_sham.times, mean(epoch_sham.data,3))
    ylabel('Fixation-related ERP')
    xlabel('Time after fixation [ms]')
    title([pID ': Sham' ], 'Interpreter','none')
    legend({EEG_epoched.chanlocs.labels})
    set(h2,'Position',[100 100 1000 660])
    linkaxes(ax)
    saveas(h2, fullfile(dir_pre, 'FRPs_eyeCA', [pID, '_fixations.png']))
    
    %% Save data
    pop_saveset(EEG_epoched, 'filename',pID, 'filepath',fullfile(dir_pre, 'FRPs_eyeCA'), 'savemode','onefile')
    
    toc
end
