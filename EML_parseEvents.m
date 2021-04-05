
% Parse triggers and other event markers for EyeMindLink EEG
%
% EEG
% --can be streamed (filename EML1_xxx) &/ recorded on SD card (LAxxxxx)
% --hardware triggers sent through the liveamp: single bit of information,
% saved in both streamed ("BV") and SD card files
% --LSL EEG & triggers direct from stimulus PC recorded in LabRececorder
% .xdf file. Note that LSL cannot accoutn for hardware delay/bluetooth delay, so although
% it would be convenient just to use the .xdf we only use it as backup.

% Trigger alignment
% --we preferentially use the hardware triggers
% --because of bluetooth connectivity dropouts, we prefer to use the SD
% card recording if it is (i) available and (ii) has sufficient triggers
% --if SD card recording is unavailable or deemed problematic (e.g. too few
% triggers) then we use the streamed recording (referred to as "BV" for
% BrainVision Recorder
% --if neither of the above are valid, we can use the XDF stream for the
% triggers

% This script synchronises hardware triggers to entries in the Trial log file
% Can deal with the following issues:
% --Some hardware triggers are missing
% --Some excess hardware triggers
% --No or so few hardware triggers that alignment cannot proceed: uses XDF
% triggers
%
%
% outputs table with events from log, datetime from log & EEG sample number
% for both streamed and SD recording
%
% EYETRACKER DATA
% read eyetracker data to find equivalent eyetracker timestamp for each
% event
% read behavioural data file to get pageread/questionpage durations (it is not simply the
% differece between successive events)
% get trial metadata such as page number, text, qType
% get correct/incorrect & MW labels
%
% **Note**: if you are reusing this code for another project (sorry) then you
% need to check your EEG sampling rate and your eyetracker sample rate.
% Ours are both at 1000Hz so I didn't need to resample.
%
% Rosy Southwell Feb 2021

%%
clear all; close all

datapath = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EyeMindLink\Data';
% sublist = % type 1: one-to-one log to trigger mapping, no hacking necessary
% events present
sublist = [19:108]; % subject 19 is the first with EEG.
beh_data = readtable('..\..\Data\EML1_allResponsesMain.csv');
beh_data = renamevars(beh_data, 'identifier', 'EVENT');
beh_data = removevars(beh_data, 'Var1');
trigSources = readtable('triggerSources.csv'); % which trigger sources are available and valid (as manually decided by Rosy based on trigger alignment in prior runs of this script)
% temp debug
sublist = [73 76 78 87 102 99]; % had NaNs for page numbers

for s = 1:length(sublist)
    pID = ['EML1_',sprintf('%03d',sublist(s))];
    fprintf('\n\n\n')
    no_sd = 0; no_bv = 0; % flags for if reading from the SD and stream fail
    if ~exist(fullfile(datapath, pID, 'EEG'),'dir')
        disp([pID ': no EEG recorded for this participant.'])
        continue
    end
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
        % logtrig.MSG = join(join(join(cellstr(string(logtrig.datetime)),logtrig.Var3,' '),join(logtrig.Var3,logtrig.Var4,' '),' '),logtrig.Var5,'  ');
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
    % initialise
    eeg_start_pc1abs = [];eeg_start_pc2abs=[];
    logtrig.eeg_sample = NaN(height(logtrig),1);
    logtrig.EEG_lag_log = NaN(height(logtrig),1);
    eeg_hardtrig=[];
    if isempty(BVfilename)
        disp('No streamed EEG file (EML1_xxx.vhdr) found. Check filenames.')
        no_bv = 1;
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
            no_bv = 1;
            
        elseif height(eeg_hardtrig)<25 % this is fairly arbitrary cutoff
            disp(['Insufficient (' num2str(height(eeg_hardtrig)) ') triggers found in streamed .vmrk!'])
            no_bv = 1;
            
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
            [x_keep, y_keep, offset, lag,hh]=align_events_diff(eeg_hardtrig.sample_corrected,milliseconds(logtrig.datetime-min(logtrig.datetime)),1000);
            figure(hh);
            title(['Hardware trigger for streamed recording (top) aligned to log (bottom): ' pID])
            saveas(hh,fullfile(datapath, pID, 'EEG','streamed_triggers.png'))
            
            
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
                no_bv = 1;
                
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
    % initialise
    logtrig.SD_lag_BV = NaN(height(logtrig),1);
    logtrig.eegSD_sample = NaN(height(logtrig),1);
    eegSD_start_pc1abs = []; eegSD_start_pc2abs = [];
    eegSD_hardtrig=[];
    if isempty(SDfilename)
        disp('No SD card file (LAxxxxx.vhdr) found. Check filenames')
        no_sd=1;
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
            no_sd=1;
        elseif height(eegSD_hardtrig)<25 % this is fairly arbitrary cutoff
            disp(['Insufficient (' num2str(height(eegSD_hardtrig)) ') triggers found in SD .vmrk!'])
            no_sd=1;
        else
            %%  align SD card hardware triggers to events in the log
            [x_keepSD, y_keepSD, offsetSD, lagSD,hh]=align_events_diff(eegSD_hardtrig.sample,milliseconds(logtrig.datetime-min(logtrig.datetime)),1000);
            figure(hh);
            title(['Hardware trigger for SDcard recording (top) aligned to log (bottom): ' pID])
            saveas(hh,fullfile(datapath, pID, 'EEG','SDcard_triggers.png'))
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
                no_sd=1;
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
    
    %%%% XDF %%%%
    %% Read XDF if neither BV or SD file is usable
    % To get events in terms of EEG sample number, we only have EEG start
    % time (as YYMMDDHHMMSSuuuuuu) and LSL triggres in terms of
    % time-since-boot of each PC. But the log has timestamps in yymmdd so
    
    %     if isempty(eeg_hardtrig) && isempty(eegSD_hardtrig)
    if no_bv && no_sd
        disp('No (usable) hardware triggers')
        if isempty(eeg_start_pc2abs) && isempty(eegSD_start_pc2abs)
            disp('Can''t find a start time for streamed or SD card EEG file. This EEG data is unusable - skipping.')
            diary off
            continue
        else
            disp('Attempting rough timing from LSL triggers')
        end
        try
            % load XDF files in both time domains
            xdf_unsync = load_xdf(fullfile(datapath, pID, [pID,'.xdf']),'HandleJitterRemoval' , false,'Verbose',false,'HandleClockSynchronization',false);
            xdf_sync = load_xdf(fullfile(datapath, pID, [pID,'.xdf']),'HandleJitterRemoval' ,false,'Verbose',false,'HandleClockSynchronization',true);
        catch
            disp('load_xdf failed - perhaps the xdf file is corrupted. ')
            disp('This EEG data is unusable - skipping.')
            diary off
            continue
        end
        % select the trigger stream
        trig_ix = find(strcmp(cellfun(@(sas) sas.info.name, xdf_unsync,'uni',false),{'eyeLink_trigger'}));
        if isempty(trig_ix)
            disp('No LSL stream named eyeLink_trigger. Will skip participant. ')
            disp('Maybe in a future iteration of this script someone may have implemented using NIRStar trig.')
            continue
        end
        xdf_trig_pc1 = xdf_sync{trig_ix};
        xdf_trig_pc2 = xdf_unsync{trig_ix};
        
        eeg_ix = find(contains(cellfun(@(sas) sas.info.name, xdf_sync,'uni',false),{'EEG'}));
        eeg_lsl_exists = ~isempty(eeg_ix);
        %         if eeg_lsl_exists
        %             % use lsl eeg
        %         else
        %% get the clock offset from the xdf
        % this is offset between boot times not between absolute times.
        try
            all_offset_tpc1 = str2double({cell2mat(xdf_trig_pc1.info.clock_offsets.offset).time})';
            all_offset_val = str2double({cell2mat(xdf_trig_pc1.info.clock_offsets.offset).value})'; % note this is pc2 time minus pc1
        catch
            disp('No clock offset recorded - alignment of EEG to log is not currently possible.')
            disp('If EEG timeseries was recorded over LSL you could use that but it will be unreliable.')
            
            disp('This EEG data is unusable - skipping.')
            diary off
            continue
        end
        mean_offset = mean(all_offset_val); % in sec
        range_offset = range(all_offset_val); % in sec
        resid_offset = all_offset_val - mean_offset; % in sec, will use later to fine tune timings
        figure(99); clf
        scatter(all_offset_tpc1,all_offset_val,'.')
        xlabel('time (s since PC2 boot)')
        ylabel('offset (boot time PC2 - PC1)')
        title([pID ' LSL clock offsets'], 'interpreter','none')
        
        %% COmpare log and LSL timestamps
        figure(22);clf
        plot(logtrig.datetime(1:min([length(xdf_trig_pc1.time_stamps) height(logtrig)])),...
            xdf_trig_pc2.time_stamps(1:min([length(xdf_trig_pc1.time_stamps) height(logtrig)])),'o' )
        title([pID ' triallog vs LSL triggers'], 'interpreter','none')
        xlabel('log timestamp')
        ylabel('LSL timestamp')
        lsl_lag_log_range = milliseconds(range(logtrig.datetime(1:min([length(xdf_trig_pc1.time_stamps) height(logtrig)])) - ...
            seconds(xdf_trig_pc2.time_stamps(1:min([length(xdf_trig_pc1.time_stamps) height(logtrig)])))'));
        %% Estimate PC1 absolute boot time
        % subtract the LSL triggers uncorrected timestamps from the log
        % timestamps, should all be the same or v similar
        % note there are some weird delays introduced later in the session
        % where lsl timestamps lag behind the log, so only take events upt o
        % event_cutoff
        event_cutoff = 25; % use earliest n events (later ones may have mroe drift)
        pc1_boot_ests = logtrig.datetime(1:event_cutoff) - ...
            seconds(xdf_trig_pc2.time_stamps(1:event_cutoff))';
        
        % how do these estimates look over time?
        figure(303);clf
        scatter(xdf_trig_pc1.time_stamps(1:event_cutoff),pc1_boot_ests(1:event_cutoff))
        title([pID ' estimation of boot time PC1'], 'interpreter','none')
        xlabel('rel timestamp (PC1)'); ylabel('estimate boot time')
        
        %  Any variability in these estimates indicates the logging and LSL
        %  timestamp functions are not getting the same clocktime as one another.
        %
        %  which estimate should we use for PC1 boot time?
        %  -- I'd say the one derived from the earliest trigger event
        %  we could actually extrapolate this back in time to timestamp (PC1) ==0 if it looks like a
        %  simple function, or take an average over several early events
        
        % % 1. earliest estimate
        % pc1_boot_abs = pc1_boot_ests(1);
        % 2. average over earliest n estimates
        pc1_boot_abs = mean(pc1_boot_ests(1:event_cutoff));
        % % 3. linear fit (check that the fn looks linear... you ma ywish to choose event_cutoff wisely here)
        %         pc1_boot_fn = polyfit(xdf_trig_pc1.time_stamps(1:event_cutoff),seconds(pc1_boot_ests-min(pc1_boot_ests)),1);
        %         figure(303);hold on;
        %         x=linspace(0,max(xdf_trig_pc1.time_stamps),100);
        %         plot(x,seconds(x*pc1_boot_fn(1)+pc1_boot_fn(2))+min(pc1_boot_ests))
        %         pc1_boot_abs = pc1_boot_fn(2)+min(pc1_boot_ests);
        
        
        %% Estimate PC2 absolute boot time
        % we know the offset between pc1 boot and pc2 boot for each sync - this is
        % the clock_offset!
        % This does seem light a tight linear relationship so let's extrapolate
        % this back in time to pc2 rel ==0
        
        % 1. get lags in terms of PC2 time
        all_offset_tpc2 = all_offset_tpc1+all_offset_val;
        % 2. fit linear func
        offset_fn_tpc2 = polyfit(all_offset_tpc2, all_offset_tpc1,1);
        
        x=linspace(0,max(all_offset_tpc2),5);
        figure(69);clf; scatter(all_offset_tpc2, all_offset_tpc1);
        hold on; plot(x,x*offset_fn_tpc2(1)+offset_fn_tpc2(2))
        % 3. extrapolate back to time (PC2) ==0
        pc2_boot_abs = pc1_boot_abs + seconds(offset_fn_tpc2(2));
        title([pID ' extrapolate back to boot time PC2'],'interpreter','none')
        xlabel('timestamp (PC2)'); ylabel('timestamp (PC1)')
        
        %% Get EEG start time in other formats. Add this to the 0-based indexing of EEG timestamps to get it in that format
        if ~isempty(eeg_start_pc2abs)
            
            eeg_start_pc2rel = seconds(eeg_start_pc2abs-pc2_boot_abs);
            eeg_start_pc1rel = offset_fn_tpc2(1)*eeg_start_pc2rel + offset_fn_tpc2(2);
            eeg_start_pc1abs = seconds(eeg_start_pc1rel)+pc1_boot_abs;
            eeg_start_pc1abs.Format=('yyyy-MM-dd HH:mm:ss.SSSSSS');
        end
        if ~isempty(eegSD_start_pc2abs)
            eegSD_start_pc2rel = seconds(eegSD_start_pc2abs-pc2_boot_abs);
            eegSD_start_pc1rel = offset_fn_tpc2(1)*eegSD_start_pc2rel + offset_fn_tpc2(2);
            eegSD_start_pc1abs = seconds(eegSD_start_pc1rel)+pc1_boot_abs;
            eegSD_start_pc1abs.Format=('yyyy-MM-dd HH:mm:ss.SSSSSS');
        end
        %% Estimate EEG sample for each logged event
        
        
        % LSL alignment disgnostics
        timinfo(s).lsl_lag_log_range = lsl_lag_log_range;
        %         if lsl_lag_log_range > 1000
        %             disp(['WARNING: LSL and log have inconsistent timing with a jitter of ' num2str(lsl_lag_log_range) 'ms'])
        %             disp('You should discard this subject from analysis.')
        %             diary off
        %
        %             continue
        %         else
        if ~isempty(eeg_start_pc2abs)
            logtrig.eeg_sample_est = round(milliseconds(logtrig.datetime - eeg_start_pc1abs));
            disp('estimated EEG-streamed sample number via LSL hacking')
        end
        if ~isempty(eegSD_start_pc2abs)
            logtrig.eegSD_sample_est =  round(milliseconds(logtrig.datetime - eegSD_start_pc1abs));
            disp('estimated EEG-SD sample number via LSL hacking')
            
        end
        %         end
        %         end
    end
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
    
    %% Get eyetracker timestamps
    msgfiles = dir(fullfile(datapath,pID,'Unpacked','*Message.csv'));
    messages=[];
    if isempty(msgfiles)
        disp(['No eyetracker Message files found! Skipping eyetracker timestamps for ' pID])
    else
        for i = 1:length(msgfiles)
            messages =  [messages;...
                readtable(fullfile(msgfiles(i).folder,msgfiles(i).name))];
        end
        
        % get eyetracker timestamp rel to eeg timestamp by matching log to ET messages
        messages_trialonsets = messages(contains(messages.text,'TRIALID ') & ~contains(messages.text,{'DriftCorrect','Recal'}),:);
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
    end
    %% Get page times and behavioural conditions
    thisbeh = beh_data(string(beh_data.participantID)==pID,:);
    thisbeh = removevars(thisbeh,'participantID');
    thisbeh = thisbeh(string(thisbeh.EVENT)~='DriftCorrect',:); % remove drift corrects / recals to avoid errors
    thisbeh = thisbeh(string(thisbeh.EVENT)~='Recal',:); % remove drift corrects / recals to avoid errors
    
    logtrig = sortrows(outerjoin(logtrig, thisbeh, 'Keys','EVENT','Type','left'));
    logtrig=renamevars(logtrig, 'EVENT_logtrig','EVENT' );
    logtrig=removevars(logtrig, 'EVENT_thisbeh');
    %% write updated event log
    writetable(logtrig,fullfile(datapath, pID, 'EEG','events.csv'))
    
    
end