% Parse triggers for outcome type 1, 2 an 3:

% EEG recorded in BV recorder &/ SD card of the LiveAmp
% Hardware triggers sent through the liveamp
% Can deal with the following issues:
% --Some hardware triggers are missing
% --Some excess hardware triggers
%
% synchronisation:
% compute timing diff between successive triggers and use timeseries
% alignment to get a matching from log to hardware triggers
% use this to identify log events with no hardware trigger and estimate EEG
% sample from the log.
% output table with events from log, datetime from log & EEG sample number
% for both streamed and SD recording

%%
clear all; close all

datapath = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EyeMindLink\Data';
% sublist = % type 1: one-to-one log to trigger mapping, no hacking necessary
% events present
sublist = [97:99];

for s = 1:length(sublist)
    pID = ['EML1_',sprintf('%03d',sublist(s))];
    fprintf('\n\n\n')
    % set up file to record some trigger diagnostics
    dfile=fullfile(datapath, pID, 'EEG',[ pID '_trigger_diagnostics.txt']);
    if exist(dfile, 'file') ; delete(dfile); end
    diary(dfile)
    diary on
    disp(['______Trigger diagnostics for ' pID '______'])
    clear logtrig eeg_hardtrig eegSD_hardtrig eeg_start eeg_start_pc1abs eeg_start_pc2abs
    
    %% get the absolute timestamps from the log file on PC1
    % earlier participants (software v3) have different log
    if sublist(s) <27
        %TODO fix this
        logtrig = readtable(fullfile(datapath,pID,[pID '_Trials.txt']) );
        % combine date and time to a datetime obj
        logtrig.Var1.Format = 'yyyy-MM-dd HH:mm:ss.SSSSSS';
        logtrig.datetime = logtrig.Var1 + logtrig.Var2;
        logtrig.Var1.Format = 'yyyy-MM-dd';
        logtrig.Properties.VariableNames{6} = 'EVENT';
        logtrig.Properties.VariableNames{7} = 'VAL';
        logtrig.MSG = join(join(join(cellstr(string(logtrig.datetime)),logtrig.Var3,' '),join(logtrig.Var3,logtrig.Var4,' '),' '),logtrig.Var5,'  ');
    else
        logtrig = readtable(fullfile(datapath,pID,[pID '_Trials.txt']) ,'Delimiter','\t');
        logtrig.Properties.VariableNames{2} = 'EVENT';
        logtrig.Properties.VariableNames{3} = 'VAL';
        % split the first column again by space delimiter
        cols = split(logtrig.Var1);
        logtrig.datetime =datetime(join( cols(:,1:2),' '));
        logtrig.datetime.Format = 'yyyy-MM-dd HH:mm:ss.SSSSSS';
        logtrig.Properties.VariableNames{1} = 'MSG';
    end
    
    % identify events during EEG recording
    on_ix = ~contains(logtrig.EVENT,'Y_');
    logtrig = logtrig(on_ix,:);
    
    %% Remove events too close for reliable trigger resolution (<5ms)
    % - take the first of the two. Of couse the first trigger is given a diff of 0 so we
    % mustn't discard that
    logtrig.diff_since_last = [0; milliseconds(diff(logtrig.datetime))];
    logtrig = logtrig(logtrig.diff_since_last(2:end)>5,:);
    
    %%%% STREAMED FILE %%%%
    %% Read EEG markers from streamed file
    BVfilename = dir(fullfile(datapath, pID, 'EEG','EML*.vmrk'));
    if isempty(BVfilename)
        disp('No streamed EEG file (EML1_xxx.vhdr) found. Check filenames')
        eeg_start_pc1abs = [];eeg_start_pc2abs=[];
        logtrig.eeg_sample = NaN(height(logtrig),1);
        logtrig.EEG_lag_log = NaN(height(logtrig),1);
    else
        temp = fileread(fullfile(datapath, pID, 'EEG', BVfilename.name));
        eeg_start = regexp(temp, 'Mk1=New Segment,,1,1,0,(?<tst>[0-9]{20})','names').tst;
        eeg_start = [eeg_start(1:4) '-' eeg_start(5:6) '-' eeg_start(7:8) ' ' eeg_start(9:10) ':' eeg_start(11:12) ':' eeg_start(13:14) '.' eeg_start(15:end)];
        eeg_start_pc2abs = datetime(eeg_start, 'Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
        
        % read in trigger file - as streamed to BV recorder
        vmrk =  readtable(fullfile(datapath, pID, 'EEG', BVfilename.name),'filetype','text','HeaderLines',11,...
            'ReadVariableNames',false,'Delimiter','comma');
        vmrk.Properties.VariableNames = {'number','comment','sample','size','channel','date'};
        eeg_hardtrig = vmrk(contains(string(vmrk.comment),"M  1"),:);
        if isempty(eeg_hardtrig)
            disp(['No triggers found in streamed .vmrk!'])
            logtrig.eeg_sample = NaN(height(logtrig),1);
            logtrig.EEG_lag_log = NaN(height(logtrig),1);
        else
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
                disp(['EEG streamed over Bluetooth: ' num2str(sum(dropped.n_dropped)) ' samples were dropped in total over ' num2str(height(dropped)) ' dropout events'])
            end
            
            
            %%  align hardware triggers to events in the log
            [x_keep, y_keep, offset, lag]=align_events_diff(eeg_hardtrig.sample_corrected,milliseconds(logtrig.datetime-min(logtrig.datetime)),1000);
            
            % store aligned events
            logtrig.eeg_sample = NaN(height(logtrig),1);
            logtrig.eeg_sample(y_keep) =  eeg_hardtrig.sample(x_keep);
            logtrig.EEG_lag_log=NaN(height(logtrig),1);
            logtrig.EEG_lag_log(y_keep) = seconds(eeg_hardtrig.PC2datetime(x_keep)-logtrig.datetime(y_keep));
            eeg_hardtrig.PC1datetime = NaT(height(eeg_hardtrig),1);
            eeg_hardtrig.PC1datetime(x_keep) = logtrig.datetime(y_keep);% get time from log
            eeg_hardtrig.val =NaN(height(eeg_hardtrig),1);
            eeg_hardtrig.val(x_keep) = logtrig.VAL(y_keep);
            eeg_hardtrig.eventLabel =cell(height(eeg_hardtrig),1);
            eeg_hardtrig.eventLabel(x_keep) = logtrig.EVENT(y_keep);
            eeg_hardtrig.EEG_lag_log = seconds(eeg_hardtrig.PC2datetime-eeg_hardtrig.PC1datetime);
            
            % estimate eeg_start_pc1abs
            eeg_start_pc1abs = mean(logtrig.datetime(y_keep)- milliseconds(eeg_hardtrig.sample(x_keep)));
            
            
            disp([num2str(height(logtrig)) ' events in log file'])
            disp([num2str(height(eeg_hardtrig)) ' triggers in the streamed .vmrk file'])
            disp(['    ' num2str(length(y_keep)) ' triggers matched to log events'])
            %% Deal w missing hardware triggers in streamed file
            logtrig.eeg_sample_est = logtrig.eeg_sample; % by default est = actual
            
            
            
            if length(y_keep) / height(logtrig) <.25
                disp(['***WARNING!!! Too few hardware triggers from stream matched to logged events.'])
                disp('Trigger repair would be unreliable, so was not attempted. Please use .xdf or SD card recording instead.')
            elseif sum(isnan(logtrig.eeg_sample))==0
            else
                % Use matched hardware triggers to find temporal alignment to log.
                % Use relative timings in log to estimate EEG sample number
                disp(['Attempting to repair ' ...
                    num2str(sum(isnan(logtrig.eeg_sample))) ' missing triggers in streamed recording using log.'])
                to_fix = find(isnan(logtrig.eeg_sample));
                good = find(~isnan(logtrig.eeg_sample));
                for i= to_fix'
                    % nearest 'good' trigger?
                    [mn ix]=min(abs(i-good));
                    % filling in backwards or forwards?
                    if i-good(ix)<0 % backwards
                        logtrig.eeg_sample_est(i) = logtrig.eeg_sample(good(ix))-sum(logtrig.diff_since_last((1+i):good(ix)));
                    end
                    if i-good(ix)>0 % forwards
                        logtrig.eeg_sample_est(i) = logtrig.eeg_sample(good(ix))+sum(logtrig.diff_since_last((good(ix)+1):i));
                    end
                end
                logtrig.eeg_sample_est(logtrig.eeg_sample_est<0) = NaN; % if the esimated timing is before EEG recording started
                if sum(isnan(logtrig.eeg_sample_est))>0
                    disp(['    Couldn''t repair all missing streamed triggers. Still missing ' num2str(sum(isnan(logtrig.eeg_sample_est)))])
                else
                    disp('    Repaired all missing streamed triggers.')
                end
            end
            
        end
    end
    
    
    
    %%%% SD CARD %%%%
    %% Read EEG markers from SD card
    % get eeg data start time in absolute time, taken from PC2 clock
    SDfilename = dir(fullfile(datapath, pID, 'EEG','LA*.vmrk'));
    if isempty(SDfilename)
        disp('No SD card file (LAxxxxx.vhdr) found. Check filenames')
        logtrig.SD_lag_BV = NaN(height(logtrig),1);
        logtrig.eegSD_sample = NaN(height(logtrig),1);
        
    else
        temp = fileread(fullfile(datapath, pID, 'EEG', SDfilename.name));
        eegSD_start = regexp(temp, 'Mk1=New Segment,,1,1,0,(?<tst>[0-9]{20})','names').tst;
        eegSD_start = [eegSD_start(1:4) '-' eegSD_start(5:6) '-' eegSD_start(7:8) ' ' eegSD_start(9:10) ':' eegSD_start(11:12) ':' eegSD_start(13:14) '.' eegSD_start(15:end)];
        eegSD_start_pc2abs = datetime(eegSD_start, 'Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
        
        % read in trigger file - as recorded onboard the LiveAmp SD card
        eegSD_hardtrig =  readtable(fullfile(datapath, pID, 'EEG', SDfilename.name),'filetype','text','HeaderLines',11,...
            'ReadVariableNames',false,'Delimiter','comma');
        
        eegSD_hardtrig.Properties.VariableNames = {'number','comment','sample','size','channel','date'};
        eegSD_hardtrig = eegSD_hardtrig(contains(string(eegSD_hardtrig.comment),"M  1"),:);
        if isempty(eegSD_hardtrig)
            disp(['No triggers found in SD card .vmrk!'])
            y_keepSD =0;
            logtrig.SD_lag_BV = NaN(height(logtrig),1);
            
        else
            %%  align SD card hardware triggers to events in the log
            [x_keepSD, y_keepSD, offsetSD, lagSD]=align_events_diff(eegSD_hardtrig.sample,milliseconds(logtrig.datetime-min(logtrig.datetime)),1000);
            
            % store aligned events
            logtrig.eegSD_sample = NaN(height(logtrig),1);
            logtrig.eegSD_sample(y_keepSD) =  eegSD_hardtrig.sample(x_keepSD);
            eegSD_hardtrig.PC1datetime = NaT(height(eegSD_hardtrig),1);
            eegSD_hardtrig.PC1datetime(x_keepSD) = logtrig.datetime(y_keepSD);    % get time from log
            eegSD_hardtrig.val =NaN(height(eegSD_hardtrig),1);
            eegSD_hardtrig.val(x_keepSD) = logtrig.VAL(y_keepSD);
            eegSD_hardtrig.eventLabel =cell(height(eegSD_hardtrig),1);
            eegSD_hardtrig.eventLabel(x_keepSD) = logtrig.EVENT(y_keepSD);
            eegSD_hardtrig.PC2datetime = milliseconds(eegSD_hardtrig.sample) +  eegSD_start_pc2abs; % this is objective from the recording
            eegSD_hardtrig.diff_since_last = [0; diff(eegSD_hardtrig.sample)];
            eegSD_hardtrig.EEG_lag_log =  seconds(eegSD_hardtrig.PC2datetime-eegSD_hardtrig.PC1datetime);
            
            
            % lag between PC2 timestamps (i.e. based on the eeg file start time)
            % for SD card rel to BV streamed data
            try
                logtrig.SD_lag_BV = seconds(milliseconds(logtrig.eeg_sample - logtrig.eegSD_sample)-(eegSD_start_pc2abs-eeg_start_pc2abs));
            end
            
            
            disp([num2str(height(eegSD_hardtrig)) ' triggers in the SD card .vmrk file'])
            disp(['    ' num2str(length(y_keepSD)) ' triggers matched to log events'])
            
            %% Deal w missing hardware triggers for SD card
            logtrig.eegSD_sample_est = logtrig.eegSD_sample;
            if length(y_keepSD) / height(logtrig) <.25
                disp(['***WARNING!!! Too few hardware triggers from SD card matched to logged events'])
                disp('Trigger repair would be unreliable, so was not attempted. Please use .xdf or streamed recording instead.')
            elseif sum(isnan(logtrig.eegSD_sample))==0
            else
                % Use matched hardware triggers to find temporal alignment to log.
                % Use relative timings in log to estimate EEG sample number
                disp(['Attempting to repair '  ...
                    num2str(sum(isnan(logtrig.eegSD_sample))) ' missing triggers in SD recording using log.'])
                to_fix = find(isnan(logtrig.eegSD_sample));
                good = find(~isnan(logtrig.eegSD_sample));
                for i= to_fix'
                    % nearest 'good' trigger?
                    [mn ix]=min(abs(i-good));
                    % filling in backwards or forwards?
                    if i-good(ix)<0 % backwards
                        logtrig.eegSD_sample_est(i) = logtrig.eegSD_sample(good(ix))-sum(logtrig.diff_since_last((1+i):good(ix)));
                    end
                    if i-good(ix)>0 % forwards
                        logtrig.eegSD_sample_est(i) = logtrig.eegSD_sample(good(ix))+sum(logtrig.diff_since_last((good(ix)+1):i));
                    end
                end
                logtrig.eegSD_sample_est(logtrig.eegSD_sample_est<0) = NaN; % if the esimated timing is before EEG recording started
                if sum(isnan(logtrig.eegSD_sample_est))>0
                    disp(['    Couldn''t repair all missing triggers on SD. Still missing ' num2str(sum(isnan(logtrig.eegSD_sample_est)))])
                else
                    disp('    Repaired all missing triggers from SD recording.')
                end
            end
        end
    end
    %
    
    %% store diagnistics
    timinfo(s).eeg_start_pc1abs = eeg_start_pc1abs;
    timinfo(s).eeg_start_pc2abs = eeg_start_pc2abs;
    timinfo(s).eeg_start_discrepancy = seconds(eeg_start_pc2abs - eeg_start_pc1abs);
    timinfo(s).SD_lag_BV_mean = nanmean(logtrig.SD_lag_BV);
    timinfo(s).SD_lag_BV_range = range(logtrig.SD_lag_BV); % there will be a lot of jitter if there are missed samples in the streamed file!
    timinfo(s).EEG_lag_log =nanmean(logtrig.EEG_lag_log );% there will be a lot of jitter if there are missed samples in the streamed file!
    timinfo(s).EEG_lag_log_range =range(logtrig.EEG_lag_log );% there will be a lot of jitter if there are missed samples in the streamed file!
    
    disp(newline)
    disp(timinfo(s))
    diary off
    %% write updated event log
    writetable(logtrig,fullfile(datapath, pID, 'EEG','events.csv'))
    
    
end