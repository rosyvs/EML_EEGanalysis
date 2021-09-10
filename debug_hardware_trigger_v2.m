% Load EEG, vmrk,  trial logging and XDf data for a quick check for EML1

% % Recordings made on 9/30
% filename, PCwith BVrec, PC with LabRec, ISI, Software stop/timer
% trigtest1005  1   1   100 T
% trigtest1006  1   1   100 S
% trigtest1007  1   -   100 T
% trigtest1008  1   2   100 s
% trigtest1009  2   1   100 T
% trigtest1010  2   1   100 S
% trigtest1011  2   2   100 T
% trigtest1012  2   2   100 S
% trigtest1013  2   2   2000 T
clear all; close all
datapath = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EyeMindLink\ScratchSpace\EEGdebug\';


sublist =[1005:1013];
xdf_invalid = {'trigtest1007','trigtest1008','trigtest1009'};
for s = 1:length(sublist)
    filename = 'trigtest';
    pID = [filename num2str(sublist(s))];
    % for trigger timing to nto get messed up, we need to disable jitter
    % removal on xdf load, otherwise it makes itt a regular stream.
    if ~any(strcmp(xdf_invalid,pID)) % yeah i forgot to record in labrecorder fro this one oops
        xdf_uncorr = load_xdf(fullfile(datapath, pID, [pID,'.xdf']),'HandleJitterRemoval' , false,'Verbose',true,'HandleClockSynchronization',false);
        xdf_corr = load_xdf(fullfile(datapath, pID, [pID,'.xdf']),'HandleJitterRemoval' ,false,'Verbose',true,'HandleClockSynchronization',true);
        
        % detect which stream is the trigger
        trig_ix = find(contains(cellfun(@(sas) sas.info.name, xdf_uncorr,'uni',false),{'RDA'}));
        xdf_trig_pc1 = xdf_uncorr{trig_ix};
        xdf_trig_pc2 = xdf_corr{trig_ix}; % the remotely-stamped timestamps are realigned to the local (PC2) clock (Sec since boot)
        % THESE SHOULD BE THE SAME for recordings made with labrecorder on PC1
        
        xdf_trig_pc1.diff_since_last = [0; diff(xdf_trig_pc1.time_stamps)'];
        xdf_trig_pc2.diff_since_last = [0; diff(xdf_trig_pc2.time_stamps)'];
    end
    %% Read EEG files
    % ___from BrainVision Recorder____
    % get eeg data start time in absolute time, taken from PC2 clock
    temp = fileread(fullfile(datapath, pID,  [pID, '.vmrk']));
    eeg_start = regexp(temp, 'Mk1=New Segment,,1,1,0,(?<tst>[0-9]{20})','names').tst;
    eeg_start = [eeg_start(1:4) '-' eeg_start(5:6) '-' eeg_start(7:8) ' ' eeg_start(9:10) ':' eeg_start(11:12) ':' eeg_start(13:14) '.' eeg_start(15:end)];
    eeg_start_pc2abs = datetime(eeg_start, 'Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
    
    % read in trigger file - as streamed to BV recorder
    eeg_hardtrig =  readtable(fullfile(datapath, pID, [pID, '.vmrk']),'filetype','text','HeaderLines',11); % TODO verify that it is always 11 lines
    eeg_hardtrig = eeg_hardtrig(contains(eeg_hardtrig.Var2,"M  1"),:);
    
    % ____from SD card____
    % get eeg data start time in absolute time, taken from PC2 clock
    % SDfilename = dir(fullfile(datapath, pID,'LA*.vmrk'));
    % temp = fileread(fullfile(datapath, pID,  SDfilename.name));
    % eegSD_start = regexp(temp, 'Mk1=New Segment,,1,1,0,(?<tst>[0-9]{20})','names').tst;
    % eegSD_start = [eegSD_start(1:4) '-' eegSD_start(5:6) '-' eegSD_start(7:8) ' ' eegSD_start(9:10) ':' eegSD_start(11:12) ':' eegSD_start(13:14) '.' eegSD_start(15:end)];
    % eegSD_start_pc2abs = datetime(eegSD_start, 'Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
    %
    % read in trigger file - as recorded onboard the LiveAmp SD card
    % eegSD_hardtrig =  readtable(fullfile(datapath, pID, SDfilename.name),'filetype','text','HeaderLines',11); % TODO verify that it is always 11 lines
    % eegSD_hardtrig = eegSD_hardtrig(contains(eegSD_hardtrig.Var2,"M  1"),:);
    
    
    if ~any(strcmp(xdf_invalid,pID)) % yeah i forgot to record in labrecorder fro this one oops
        
        %% get the clock offset from the xdf -- NOTE THIS DOESNT
        % this is offset between boot times not between absolute times.
        % there are multiple measurements of the offset = we need to decide how to deal w that!
        % - Simplest is to compute mean offsest, use as a constant, then consider finer adjustment later w the residuals
        % - clock_offsets are the same in teh corrected and uncorrected files -
        % its just part of the XDF header
        all_offset_tpc1 = str2double({cell2mat(xdf_trig_pc1.info.clock_offsets.offset).time})';
        all_offset_val = str2double({cell2mat(xdf_trig_pc1.info.clock_offsets.offset).value})'; % note this is pc2 time minus pc1
        
        mean_offset = mean(all_offset_val); % in sec; this is approximately the diff between boot times
        range_offset = range(all_offset_val); % in sec
        resid_offset = all_offset_val - mean_offset; % in sec, will use later to fine tune timings
        figure(99); clf
        scatter(all_offset_tpc1,all_offset_val,'.')
        xlabel('time (s since PC2 boot)')
        ylabel('offset (boot time PC2 - PC1)')
        title([pID ' LSL clock offsets'], 'interpreter','none')
    end
    %% get the absolute timestamps from the log file on PC1
    logname = dir(fullfile(datapath,pID,'*.txt'));
    logtrig = readtable(fullfile(logname.folder, logname.name) ,'Delimiter','space','HeaderLines', 1);
    logtrig = logtrig(:,[1 2 8]);
    logtrig.Properties.VariableNames{3} = 'VAL';
    % split the first column again by space delimiter
    logtrig.datetime = logtrig.Var1 + logtrig.Var2;
    logtrig.datetime.Format = 'yyyy-MM-dd HH:mm:ss.SSSSSS';
    
    
    figure(22);clf
    plot(logtrig.datetime(1:length(xdf_trig_pc1.time_stamps)),xdf_trig_pc2.time_stamps,'o' )
    title([pID ' triallog vs LSL triggers'], 'interpreter','none')
    xlabel('log timestamp')
    ylabel('LSL timestamp')
    
    if ~any(strcmp(xdf_invalid,pID)) % yeah i forgot to record in labrecorder fro this one oops
        
        %% Estimate PC1 absolute boot time from LSL clock offsets
        % subtract the LSL triggers uncorrected timestamps from the log
        % timestamps, the estimates should all be the same or v similar
        % note there are some weird delays introduced later in the session
        % where lsl timestamps lag behind the log, so only take events upt o
        % event_cutoff
        event_cutoff = 1000;
        pc1_boot_ests = logtrig.datetime(1:event_cutoff) - ...
            seconds(xdf_trig_pc1.time_stamps(1:event_cutoff))';
        
        % how do these estimates look over time?
        figure(303);clf
        scatter(xdf_trig_pc1.time_stamps(1:event_cutoff),pc1_boot_ests(1:event_cutoff))
        title([pID ' estimation of boot time PC1'], 'interpreter','none')
        xlabel('rel timestamp (PC1)'); ylabel('estimate boot time')
        
        % % 1. earliest estimate
        %pc1_boot_abs = pc1_boot_ests(1);
        
        % % 2. average over multiple estimates
        % pc1_boot_abs = mean(pc1_boot_ests(1:event_cutoff));
        
        % 3. linear fit * only appropriate if the pc1 ests vary linearly over
        % time & we believe this holds back to before the experiemnt
        %started...
        
        pc1_boot_fn = polyfit(xdf_trig_pc1.time_stamps(1:event_cutoff),datenum(pc1_boot_ests),1);
        figure(303);hold on;
        x=linspace(0,max(xdf_trig_pc1.time_stamps),5);
        plot(x,datetime(x*pc1_boot_fn(1)+pc1_boot_fn(2),'ConvertFrom','datenum'))
        pc1_boot_abs = datetime(pc1_boot_fn(2),'ConvertFrom','datenum','Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
        
        % get LSL trigs in terms of PC1 boot
        xdf_trig_pc1.PC1datetime = seconds(xdf_trig_pc1.time_stamps) + pc1_boot_abs;
        %
        % % FIND INDEX DELAY BETWEEN LSL and event log - this shold be 0, if not,
        % % then there are some extra or missing triggers in one or the other
        % LSLtrig_lag = finddelay(xdf_trig_pc1.time_series',...
        %     logtrig.VAL);
        
        %% Estimate PC2 absolute boot time
        % we know the offset between pc1 boot and pc2 boot for each sync - this is
        % the clock_offset!
        
        % 1. get lags in terms of PC2 time
        all_offset_tpc2 = all_offset_tpc1+all_offset_val;
        % 2. fit linear func
        offset_fn_tpc2 = polyfit(all_offset_tpc2, all_offset_tpc1,1);
        
        x=linspace(0,max(all_offset_tpc2),5);
        figure(69);clf; scatter(all_offset_tpc2, all_offset_tpc1);
        hold on; plot(x,x*offset_fn_tpc2(1)+offset_fn_tpc2(2))
        % 3. extrapolate back to time (PC2) ==0
        pc2_boot_abs = pc1_boot_abs + seconds(offset_fn_tpc2(2));
        title([pID ' estimation of boot time PC2'],'interpreter','none')
        xlabel('timestamp (PC2)'); ylabel('timestamp (PC1)')
    end
    
    %% Get EEG start time in other formats. Add this to the 0-based indexing of EEG timestamps to get it in that format
    eeg_start_pc2rel = seconds(eeg_start_pc2abs-pc2_boot_abs);
    eeg_start_pc1rel = offset_fn_tpc2(1)*eeg_start_pc2rel + offset_fn_tpc2(2);
    eeg_start_pc1abs = seconds(eeg_start_pc1rel)+pc1_boot_abs;
    eeg_start_pc1abs.Format=('yyyy-MM-dd HH:mm:ss.SSSSSS');
    
    % eegSD_start_pc2rel = seconds(eegSD_start_pc2abs-pc2_boot_abs);
    % eegSD_start_pc1rel = offset_fn_tpc2(1)*eegSD_start_pc2rel + offset_fn_tpc2(2);
    % eegSD_start_pc1abs = seconds(eegSD_start_pc1rel)+pc1_boot_abs;
    % eegSD_start_pc1abs.Format=('yyyy-MM-dd HH:mm:ss.SSSSSS');
    
    %% compute diffs between the logged events and hardware triggers from the BV recordign
    if ~any(strcmp(xdf_invalid,pID)) % yeah i forgot to record in labrecorder fro this one oops
        
        eeg_hardtrig.PC1datetime_lslest = milliseconds(eeg_hardtrig.Var3) +  eeg_start_pc1abs; % this estimated via LSL
    end
    eeg_hardtrig.PC2datetime = milliseconds(eeg_hardtrig.Var3) +  eeg_start_pc2abs; % this is objective from the recording
    eeg_hardtrig.diff_since_last = [0; diff(eeg_hardtrig.Var3)];
    logtrig.diff_since_last = [0; milliseconds(diff(logtrig.datetime))];
    
    % here the triggers are so close together and all come through so jus
    % tassume log[i] crorrepsonds to eeg[i]
    for i =1:height(eeg_hardtrig)
        logtrig.hardware_trig_ix(i) = i;
        
    end
    if ~any(strcmp(xdf_invalid,pID)) % yeah i forgot to record in labrecorder fro this one oops
        
        logtrig.hardware_lag_PC1est= seconds(eeg_hardtrig.PC1datetime_lslest-logtrig.datetime);
        eeg_hardtrig.hardware_lag_PC1est= seconds(eeg_hardtrig.PC1datetime_lslest-logtrig.datetime);
    end
    eeg_hardtrig.PC1datetime = logtrig.datetime;
    logtrig.hardware_lag_log =  seconds(eeg_hardtrig.PC2datetime-eeg_hardtrig.PC1datetime);
    %     eegSD_hardtrig.eventLabel = logtrig.VAL;
    if ~any(strcmp(xdf_invalid,pID)) % yeah i forgot to record in labrecorder fro this one oops
        
        % difference between log file timings and EEG timings based on the PC1
        % boot estimate
        eeg_hardtrig.PC1esterror = seconds(eeg_hardtrig.PC1datetime -  eeg_hardtrig.PC1datetime_lslest);
        % logtrig.hardware_lag_PC1est should be pretty consistent over the
        % session. The magnitude of this difference includes any error in the estimate of PC1_boot_abs.
        % BUT The jitter in this delay indicates jitter in timing of
        % hardware trigger sendding relative to trial logging; which is a
        % bigegr issue.
        
        hardware_lag_PC1est_jitter = range(logtrig.hardware_lag_PC1est); % in sec
    end
    %this is how much later the PC2 timestamps (from the hardware triggers)
    %are than PC1 log timestamps
    % This gives another estimate of PC2-PC1 clock offsets but between the absolute clocks (assuming negligible hardware delay)
    eeg_hardtrig.pc2_lag_pc1 = seconds(eeg_hardtrig.PC2datetime-eeg_hardtrig.PC1datetime);
    
    % %% compute diffs between the logged events and hardware triggers from the BV recordign
    % eegSD_hardtrig.PC1datetime_lslest = milliseconds(eegSD_hardtrig.Var3) +  eegSD_start_pc1abs; % this estimated via LSL
    % eegSD_hardtrig.PC2datetime = milliseconds(eegSD_hardtrig.Var3) +  eegSD_start_pc2abs; % this is objective from the recording
    % eegSD_hardtrig.diff_since_last = [0; diff(eegSD_hardtrig.Var3)];
    % logtrig.diff_since_last = [0; milliseconds(diff(logtrig.datetime))];
    %
    % % here the triggers are so close together and all come through so jus
    % % tassume log[i] crorrepsonds to eegSD[i]
    % for i =1:length(eegSD_hardtrig.PC1datetime_lslest)
    %     logtrig.SDhardware_trig_ix(i) = i;
    %
    % end
    %     logtrig.SDhardware_lag_PC1est= seconds(eegSD_hardtrig.PC1datetime_lslest-logtrig.datetime);
    %    eegSD_hardtrig.hardware_lag_PC1est= seconds(eegSD_hardtrig.PC1datetime_lslest-logtrig.datetime);
    %     eegSD_hardtrig.PC1datetime = logtrig.datetime;
    %     eegSDSD_hardtrig.eventLabel = logtrig.VAL;
    %
    % eegSD_hardtrig.PC1esterror = seconds(eegSD_hardtrig.PC1datetime -  eegSD_hardtrig.PC1datetime_lslest);
    % hardware_lag_PC1est_jitter = range(logtrig.hardware_lag_PC1est); % in sec
    % eegSD_hardtrig.pc2_lag_pc1 = seconds(eegSD_hardtrig.PC2datetime-eegSD_hardtrig.PC1datetime);
    % eegSD_hardtrig.SD_lag_BV = seconds(eegSD_hardtrig.PC2datetime-eeg_hardtrig.PC2datetime);
    if ~any(strcmp(xdf_invalid,pID)) % yeah i forgot to record in labrecorder fro this one oops
        
        %% compute diffs between the logged events and LSL events
        % here the triggers are so close together and all come through so jus
        % tassume log[i] crorrepsonds to eegSD[i]
        for i =1:length(xdf_trig_pc1.PC1datetime)
            
            logtrig.lsl_trig_ix(i) = i;
        end
        logtrig.lsl_log_delay= seconds(xdf_trig_pc1.PC1datetime' - logtrig.datetime);
    end
    %% Plot EEG and triggers from log file and vmrk
    % figure(3); clf;
    %
    % % plot vertical bars for the triggers in the log
    % ax1=subplot(4,1,1);
    % plot([logtrig.datetime'; logtrig.datetime'],...
    %     repmat(ylim',1,length(logtrig.datetime)))
    % %     plot(logtrig.datetime, logtrig.VAL, 's')
    % title('log file')
    % xlabel('time')
    % ylabel('trigger value')
    %
    % % plot vertical bars for the triggers in the LSL triggers
    % ax2=subplot(4,1,2);
    % plot([xdf_trig_pc1.PC1datetime; xdf_trig_pc1.PC1datetime],...
    %     repmat(ylim',1,length(xdf_trig_pc1.PC1datetime)))
    % %     plot(logtrig.datetime, logtrig.VAL, 's')
    % title('lsl')
    % xlabel('time')
    % ylabel('trigger value')
    %
    % % plot vertical bars for the thardware triggers
    % ax3=subplot(4,1,3);
    % plot([eeg_hardtrig.PC1datetime'; eeg_hardtrig.PC1datetime'],...
    %     repmat(ylim',1,length(eeg_hardtrig.PC1datetime)))
    % %     plot(logtrig.datetime, logtrig.VAL, 's')
    % title('BV hardware triggers')
    % xlabel('time')
    % ylabel('trigger value')
    
    % % plot vertical bars for the thardware triggers
    % ax4=subplot(4,1,4);
    % plot([eegSD_hardtrig.PC1datetime'; eegSD_hardtrig.PC1datetime'],...
    %     repmat(ylim',1,length(eegSD_hardtrig.PC1datetime)))
    % %     plot(logtrig.datetime, logtrig.VAL, 's')
    % title('SD card hardware triggers')
    % xlabel('time')
    % ylabel('trigger value')
    
    % linkaxes([ax1 ax2 ax3 ],'x')
    
    % sgtitle([pID ' aligned to PC1 clock'], 'interpreter','none')
    
    %% print some diagnostics
    fprintf('\n\n\n')
    disp(['______Trigger diagnistics for ' pID '______'])
    disp(['From the log file, we expected ' num2str(height(logtrig)) ' events'])
    disp(['    found ' num2str(length(xdf_trig_pc1.time_series)) ' triggers in the XDF file'])
    disp(['    found ' num2str(height(eeg_hardtrig)) ' triggers in the BV .vmrk file'])
    % disp(['    found ' num2str(height(eegSD_hardtrig)) ' triggers in the SD card .vmrk file'])
    
    %% store diagnistics
    
    if ~any(strcmp(xdf_invalid,pID)) % yeah i forgot to record in labrecorder fro this one oops
        timinfo(s).lsl_offset_jitter = range_offset;
        timinfo(s).pc1_boot_abs = pc1_boot_abs;
        timinfo(s).pc2_boot_abs = pc2_boot_abs;
        timinfo(s).pc1_boot_est_range = seconds(pc1_boot_ests(end) - pc1_boot_ests(1));
                timinfo(s).lsl_lag_log = nanmean( logtrig.lsl_log_delay);
        timinfo(s).lsl_lag_log_jitter = range( logtrig.lsl_log_delay);
        timinfo(s).eeg_start_pc1abs = eeg_start_pc1abs;
        timinfo(s).eeg_start_pc2abs = eeg_start_pc2abs;
        timinfo(s).eeg_start_discrepancy = seconds(eeg_start_pc2abs - eeg_start_pc1abs);
        timinfo(s).hardware_lag_PC1est = nanmean(logtrig.hardware_lag_PC1est);
        % timinfo(s).SDhardware_lag_PC1est = nanmean(logtrig.SDhardware_lag_PC1est);
        timinfo(s).hardware_lag_PC1est_jitter = range(logtrig.hardware_lag_PC1est);
    else
          timinfo(s).lsl_offset_jitter = NaN;
        timinfo(s).pc1_boot_abs = NaN;
        timinfo(s).pc2_boot_abs = NaN;
        timinfo(s).pc1_boot_est_range = NaN;
                        timinfo(s).lsl_lag_log = NaN;
        timinfo(s).lsl_lag_log_jitter = NaN;
        timinfo(s).eeg_start_pc1abs = NaN;
        timinfo(s).eeg_start_pc2abs = NaN;
        timinfo(s).eeg_start_discrepancy = NaN;
        timinfo(s).hardware_lag_PC1est = NaN;
        % timinfo(s).SDhardware_lag_PC1est = NaN;
        timinfo(s).hardware_lag_PC1est_jitter = NaN;
    end
    % timinfo(s).SD_lag_BV_mean = nanmean(eegSD_hardtrig.SD_lag_BV);
    % timinfo(s).SD_lag_BV_jitter = range(eegSD_hardtrig.SD_lag_BV);
    
    % timinfo(s).SDhardware_lag_PC1est_jitter =range(logtrig.SDhardware_lag_PC1est);
    timinfo(s).hardware_lag_log =nanmean(eeg_hardtrig.pc2_lag_pc1 );
    timinfo(s).hardware_lag_log_jitter =range(eeg_hardtrig.pc2_lag_pc1 );
    
    disp 'EEG trigger ISI'
    disp([num2str(min(eeg_hardtrig.diff_since_last(2:end))) ' - ' num2str(max(eeg_hardtrig.diff_since_last(2:end))) 'ms'])
    %
    % disp 'EEG SD trigger ISI'
    % disp([num2str(min(eegSD_hardtrig.diff_since_last(2:end))) ' - ' num2str(max(eegSD_hardtrig.diff_since_last(2:end))) 'ms'])
    
    disp 'actual logged ISI'
    disp([num2str(min(logtrig.diff_since_last(2:end))) ' - ' num2str(max(logtrig.diff_since_last(2:end))) 'ms'])
    
    
    disp(timinfo(s))
    clear xdf_trig_pc1 xdf_trig_pc2 logtrig eeg_hardtrig
    
end
