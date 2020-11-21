% Fix trigger issues for outcome type 2:

% EEG recorded in BV recorder &/ SD card of the LiveAmp
% Hardware triggers sent through the liveamp
% more triggers than expected came through on hardware
% (e.g. experiment restarted, cable shorts)
% synchronisation:
% compute timing diff between successive triggers and use timeseries
% alignment to get a matching from log to hardware triggers
% use this to remove superfluous hardware triggers
% output table with events from log, datetime from log & EEG sample number
% for both streamed and SD recording

%%
clear all; close all

datapath = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EyeMindLink\Data';

sublist = [56 60 63 66 67];
sublist=67
for s = 1:length(sublist)
    pID = ['EML1_',sprintf('%03d',sublist(s))];
    
    %% Read EEG markers from streamed file
    % get eeg data start time in absolute time, taken from PC2 clock
    BVfilename = dir(fullfile(datapath, pID, 'EEG','EML*.vmrk'));
    temp = fileread(fullfile(datapath, pID, 'EEG', BVfilename.name));
    eeg_start = regexp(temp, 'Mk1=New Segment,,1,1,0,(?<tst>[0-9]{20})','names').tst;
    eeg_start = [eeg_start(1:4) '-' eeg_start(5:6) '-' eeg_start(7:8) ' ' eeg_start(9:10) ':' eeg_start(11:12) ':' eeg_start(13:14) '.' eeg_start(15:end)];
    eeg_start_pc2abs = datetime(eeg_start, 'Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
    
    % read in trigger file - as streamed to BV recorder
    vmrk =  readtable(fullfile(datapath, pID, 'EEG', BVfilename.name),'filetype','text','HeaderLines',11,...
        'ReadVariableNames',false,'Delimiter','comma');
    vmrk.Properties.VariableNames = {'number','comment','sample','size','channel','date'};
    eeg_hardtrig = vmrk(contains(vmrk.comment,"M  1"),:);
    eeg_hardtrig.PC2datetime = milliseconds(eeg_hardtrig.sample) +  eeg_start_pc2abs; % this is objective from the recording
    eeg_hardtrig.diff_since_last = [0; diff(eeg_hardtrig.sample)];
    
    % read dropped samples - we will need this for aligning later~!
    % BV reords number of dropped samples, but it does not increment EEG or
    % vmrk sample number by the dropped samples. This will pose an issue
    % for aligning to the trigger log.
    dropped = vmrk(contains(vmrk.comment,"LostSamples"),:);
            eeg_hardtrig.sample_corrected = eeg_hardtrig.sample; % by defaul corrected = orig

    if ~isempty(dropped)
        bob = regexp(dropped.comment,'LostSamples: (\d+)','once','tokens');
        dropped.n_dropped = str2double([bob{:}])';
        for d=1:height(dropped)
            afters = eeg_hardtrig.sample>dropped.sample(d);
            eeg_hardtrig.sample_corrected(afters) = eeg_hardtrig.sample_corrected(afters) + dropped.n_dropped(d);
        end
    end
    
    %% Read EEG markers from SD card
    % get eeg data start time in absolute time, taken from PC2 clock
    SDfilename = dir(fullfile(datapath, pID, 'EEG','LA*.vmrk'));
    temp = fileread(fullfile(datapath, pID, 'EEG', SDfilename.name));
    eegSD_start = regexp(temp, 'Mk1=New Segment,,1,1,0,(?<tst>[0-9]{20})','names').tst;
    eegSD_start = [eegSD_start(1:4) '-' eegSD_start(5:6) '-' eegSD_start(7:8) ' ' eegSD_start(9:10) ':' eegSD_start(11:12) ':' eegSD_start(13:14) '.' eegSD_start(15:end)];
    eegSD_start_pc2abs = datetime(eegSD_start, 'Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
    
    % read in trigger file - as recorded onboard the LiveAmp SD card
    eegSD_hardtrig =  readtable(fullfile(datapath, pID, 'EEG', SDfilename.name),'filetype','text','HeaderLines',11,...
        'ReadVariableNames',false,'Delimiter','comma');
    eegSD_hardtrig.Properties.VariableNames = {'number','comment','sample','size','channel','date'};
    eegSD_hardtrig = eegSD_hardtrig(contains(eegSD_hardtrig.comment,"M  1"),:);
    eegSD_hardtrig.PC2datetime = milliseconds(eegSD_hardtrig.sample) +  eeg_start_pc2abs; % this is objective from the recording
    eegSD_hardtrig.diff_since_last = [0; diff(eegSD_hardtrig.sample)];
    
    %% get the absolute timestamps from the log file on PC1
    % earlier participants (software v3) have different log
    if sublist(s) <27
        logtrig = readtable(fullfile(datapath,pID,[pID '_Trials.txt']) );
        % combine date and time to a datetime obj
        logtrig.Var1.Format = 'yyyy-MM-dd HH:mm:ss.SSSSSS';
        logtrig.datetime = logtrig.Var1 + logtrig.Var2;
        logtrig.Var1.Format = 'yyyy-MM-dd';
        logtrig.Properties.VariableNames{7} = 'VAL';
    else
        logtrig = readtable(fullfile(datapath,pID,[pID '_Trials.txt']) ,'Delimiter','\t');
        logtrig.Properties.VariableNames{3} = 'VAL';
        logtrig.Properties.VariableNames{2} = 'MSG';
        % split the first column again by space delimiter
        cols = split(logtrig.Var1);
        logtrig.datetime =datetime(join( cols(:,1:2),' '));
        logtrig.datetime.Format = 'yyyy-MM-dd HH:mm:ss.SSSSSS';
        
    end
    
    % identify events during EEG recording
    on_ix = ~contains(logtrig.MSG,'Y_');
    logtrig = logtrig(on_ix,:);
    
        %% Remove events too close for reliable trigger resolution (<10ms) - take the
    % first of the two. Of couse the first one is given a diff of 0 so we
    % mustn't discard that
    logtrig.diff_since_last = [0; milliseconds(diff(logtrig.datetime))];
    logtrig = logtrig(logtrig.diff_since_last(2:end)>10,:);
    %%  align hardware triggers to events in the log
    [x_keep, y_keep, offset, lag]=align_events_diff(eeg_hardtrig.sample_corrected,milliseconds(logtrig.datetime-min(logtrig.datetime)),100);
    
    % store aligned events
    logtrig.eeg_sample = NaN(height(logtrig),1);
    logtrig.eeg_sample(y_keep) =  eeg_hardtrig.sample(x_keep);
    
    eeg_hardtrig.PC1datetime = NaT(height(eeg_hardtrig),1);
    eeg_hardtrig.PC1datetime(x_keep) = logtrig.datetime(y_keep);% get time from log
    eeg_hardtrig.val =NaN(height(eeg_hardtrig),1);
    eeg_hardtrig.val(x_keep) = logtrig.VAL(y_keep);
    eeg_hardtrig.eventLabel =cell(height(eeg_hardtrig),1);
    eeg_hardtrig.eventLabel(x_keep) = logtrig.MSG(y_keep);
    eeg_hardtrig.hardware_lag_log = seconds(eeg_hardtrig.PC2datetime-eeg_hardtrig.PC1datetime);
    
    % estimate eeg_start_pc1abs
    eeg_start_pc1abs = mean(logtrig.datetime(y_keep)- milliseconds(eeg_hardtrig.sample(x_keep)));
    
    %%  align SD card hardware triggers to events in the log
    [x_keepSD, y_keepSD, offsetSD, lagSD]=align_events_diff(eegSD_hardtrig.sample,milliseconds(logtrig.datetime-min(logtrig.datetime)),100);
    
    % store aligned events
    logtrig.eegSD_sample = NaN(height(logtrig),1);
    logtrig.eegSD_sample(y_keep) =  eegSD_hardtrig.sample(x_keep);
    
    eegSD_hardtrig.PC1datetime = NaT(height(eegSD_hardtrig),1);
    eegSD_hardtrig.PC1datetime(x_keep) = logtrig.datetime(y_keep);    % get time from log
    eegSD_hardtrig.val =NaN(height(eegSD_hardtrig),1);
    eegSD_hardtrig.val(x_keep) = logtrig.VAL(y_keep);
    eegSD_hardtrig.eventLabel =cell(height(eegSD_hardtrig),1);
    eegSD_hardtrig.eventLabel(x_keep) = logtrig.MSG(y_keep);
    eegSD_hardtrig.hardware_lag_log =  seconds(eegSD_hardtrig.PC2datetime-eegSD_hardtrig.PC1datetime);
    
    
    % lag between PC2 timestamps (i.e. based on the eeg file start time)
    % for SD card rel to BV streamed data
    logtrig.SD_lag_BV = seconds(milliseconds(logtrig.eeg_sample - logtrig.eegSD_sample)-(eegSD_start_pc2abs-eeg_start_pc2abs));
    
    %% print some diagnostics
    fprintf('\n\n\n')
    disp(['______Trigger diagnistics for ' pID '______'])
    disp(['From the log file, we expected ' num2str(height(logtrig)) ' events'])
    disp(['    found ' num2str(height(eeg_hardtrig)) ' triggers in the BV .vmrk file'])
    disp(['       of which ' num2str(length(y_keep)) ' matched to log events'])
    if length(y_keep)<height(logtrig)
        disp(['    WARNING: missing ' (num2str(height(logtrig))-num2str(length(y_keep))) ' triggers'])
    end
    disp(['    found ' num2str(height(eegSD_hardtrig)) ' triggers in the SD card .vmrk file'])
    disp(['       of which ' num2str(length(y_keepSD)) ' matched to log events'])
    if length(y_keepSD)<height(logtrig)
        disp(['    WARNING: missing ' (num2str(height(logtrig))-num2str(length(y_keepSD))) ' triggers'])
    end
    %% store diagnistics
    timinfo(s).eeg_start_pc1abs = eeg_start_pc1abs;
    timinfo(s).eeg_start_pc2abs = eeg_start_pc2abs;
    timinfo(s).eeg_start_discrepancy = seconds(eeg_start_pc2abs - eeg_start_pc1abs);
    timinfo(s).SD_lag_BV_mean = nanmean(logtrig.SD_lag_BV);
    timinfo(s).SD_lag_BV_jitter = range(logtrig.SD_lag_BV); % there will be a lot of jitter if there are missed samples in the streamed file!
    timinfo(s).hardware_lag_log =nanmean(eeg_hardtrig.hardware_lag_log );
    timinfo(s).hardware_lag_log_jitter =range(eeg_hardtrig.hardware_lag_log );
    
    disp(timinfo(s))
end