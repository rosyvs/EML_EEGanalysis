% fixations sham vs word

% eyeCA cleaned vs uncleaned data
% These FRPS have NOT been overlap corrected

clear all; close all
eeglab nogui % sets path defaults

%%%%%%%%%%%

do_processing = 1;
%%%%%%%%%%%


% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
sublist = [100];
exclude = [20:26, 54, 73, 77]; % Subj to exclude because no eeg or no trigger etc. 54 i fixable - check message file 6
sublist = sublist(~ismember(sublist,exclude) & ismember(sublist,find(hasTriggerList.sdcard==1)));

dir_events = '/Volumes/Blue1TB/EyeMindLink/Proessed/events/';
dir_pre = fullfile('..','..','Data','EEG_processed') ;
dir_in = fullfile(dir_pre, 'preprocessed_set');
dir_fif = fullfile(dir_pre, 'preprocessed_fif');

h_compfig=figure('name','components');
h1=figure('name','Saccade-related potentials before and after OPTICAT-ICA');
h2=figure('name','Fixation-related potentials before and after OPTICAT-ICA');
h3=figure('name','Blink-related potentials before and after OPTICAT-ICA');

analysis = 'FRP_reading_sham';
mkdir( dir_pre, analysis)
mkdir(fullfile(dir_pre, analysis, 'subavg'))
commands = {};
%%
if do_processing
    for s = 1:length(sublist)
        tic;
        close all
        pID = ['EML1_',sprintf('%03d',sublist(s))];
        
        EEG = pop_loadset([pID '.set'], dir_in);
        events = struct2table(EEG.event);
        eye_events = events(contains(events.type, 'eye'),:);
        %% info
        % Read info txt to determine whether EEG+triggers are from SD card
        % (default) or streamed (backup)
        triginfo = readtxtfile(fullfile(dir_fif, [pID '-info.txt']));
        
        % read events.csv for triggers and descriptions
        logtrig = readtable(fullfile(dir_events,[pID '_events.csv']));
        % copy the correct EEGsample column for use depending on triginfo
        if contains(triginfo, 'LA0','IgnoreCase',false)
            logtrig.eeg_use_sample = logtrig.eegSD_sample_est;
        else
            logtrig.eeg_use_sample = logtrig.eeg_sample_est;
        end
        
        %% make new event structure with fixations by reading vs sham
        task_events = logtrig(:,{'EVENT_logtrig','VAL','eeg_use_sample','duration_sec'});
        task_events = renamevars(task_events,{'EVENT_logtrig','eeg_use_sample','duration_sec'},{'type','latency','duration'});
        task_events.duration = task_events.duration * EEG.srate; % into samples
        task_events = task_events(task_events.VAL==7|task_events.VAL==20,:);
        
        task_events.type(task_events.VAL==7) = {'reading'};
        task_events.type(task_events.VAL==20) = {'sham'};
        task_events.urevent=NaN(height(task_events),1);
        task_events = removevars(task_events,'VAL');
        
        for e = 1:height(task_events)
            sel = find((eye_events.latency >= task_events.latency(e)) & (eye_events.latency < task_events.latency(e)+task_events.duration(e) ));
            eye_events.type(sel) = strcat(eye_events.type(sel), '_WITHIN_', repmat(task_events.type(e),length(sel),1));
        end
        
        eye_events = sortrows([eye_events; task_events(:,{'type','latency','duration'})],'latency');
        
        % remove NaN latency events
        eye_events = eye_events(~isnan(eye_events.latency),:);
        
        EEG.event = table2struct(eye_events);
        EEG = eeg_checkset(EEG);
        
        %% epoch into sham/reading trials
        % Preprocess data
        EEGp = pop_eegfiltnew(EEG,0.1,100);
        EEGp = pop_resample(EEGp, 200);
        EEGp = eeg_checkset(EEGp);
        [EEG_epoched, epoch_ix] = pop_epoch(EEGp, {'fixation_either_eye_WITHIN_reading','fixation_either_eye_WITHIN_sham'},[-0.2,0.8]);
        epoch_reading = pop_rmbase(pop_epoch(EEGp,{'fixation_either_eye_WITHIN_reading'},[-0.2 0.8]),[-100 0]);
        epoch_sham = pop_rmbase(pop_epoch(EEGp,{'fixation_either_eye_WITHIN_sham'},[-0.2 0.8]),[-100 0]);
        % subsample to equate n epochs
        [n_epoch_bal, which_smaller]= min([length(epoch_sham.epoch), length(epoch_reading.epoch)]);
        if which_smaller == 1
            v = 1:length(epoch_reading.epoch);
            epoch_subix = v(randperm(n_epoch_bal));
            epoch_reading = pop_select(epoch_reading, 'trial',epoch_subix);
        else
            v = 1:length(epoch_sham.epoch);
            epoch_subix = v(randperm(n_epoch_bal));
            epoch_sham = pop_select(epoch_sham, 'trial',epoch_subix);
            
        end
        
        
        %% Plot epochs
        % pop_erpimage(EEG_epoched,1)
        
        
        subavg.reading = mean(epoch_reading.data,3);
        subavg.sham = mean(epoch_sham.data,3);
        subavg.time = epoch_sham.times;
        subavg.channels = epoch_sham.chanlocs.labels;
        
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
        saveas(h2, fullfile(dir_pre, analysis, [pID '_'  analysis '.png']))
        
        %% Save data
        pop_saveset(EEG_epoched, 'filename',pID, 'filepath',fullfile(dir_pre, analysis), 'savemode','onefile')
        
        save(fullfile(dir_pre, analysis,'subavg', [pID '.mat']), 'subavg');
        %     % EEGLAB STUDY commands
        %     command_i = {'index' s 'load' fullfile(fullfile(dir_pre, analysis, [pID '.set'])) 'subject' pID };
        %     commands{s} = command_i;
        % [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG_epoched);
        toc
    end
end



% %TODO group level analyses
% [STUDY, ALLEEG] = std_editset( [], [], 'name','EML',...
%         'task', analysis,...
%         'filename', [analysis '.study'],'filepath', fullfile(dir_pre,analysis),...
%         'commands', { ...
%         { 'index' 1 'load' fullfile(filepath, 's02','syn02-s253-clean.set') 'subject' 'S02' 'condition' 'synonyms' }, ...
%         { 'index' 2 'load' fullfile(filepath, 's05', 'syn05-s253-clean.set') 'subject' 'S05' 'condition' 'synonyms' }, ...
%         { 'index' 3 'load' fullfile(filepath, 's07', 'syn07-s253-clean.set') 'subject' 'S07' 'condition' 'synonyms' }, ...
%         { 'index' 4 'load' fullfile(filepath, 's08', 'syn08-s253-clean.set') 'subject' 'S08' 'condition' 'synonyms' }, ...
%         { 'index' 5 'load' fullfile(filepath, 's10', 'syn10-s253-clean.set') 'subject' 'S10' 'condition' 'synonyms' }, ...
%         { 'index' 6 'load' fullfile(filepath, 's02', 'syn02-s254-clean.set') 'subject' 'S02' 'condition' 'non-synonyms' }, ...
%         { 'index' 7 'load' fullfile(filepath, 's05', 'syn05-s254-clean.set') 'subject' 'S05' 'condition' 'non-synonyms' }, ...
%         { 'index' 8 'load' fullfile(filepath, 's07', 'syn07-s254-clean.set') 'subject' 'S07' 'condition' 'non-synonyms' }, ...
%         { 'index' 9 'load' fullfile(filepath, 's08', 'syn08-s254-clean.set') 'subject' 'S08' 'condition' 'non-synonyms' }, ...
%         { 'index' 10 'load' fullfile(filepath, 's10', 'syn10-s254-clean.set') 'subject' 'S10' 'condition' 'non-synonyms' }, ...
%     { 'dipselect' 0.15 } });