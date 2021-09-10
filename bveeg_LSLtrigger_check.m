% Load EEG, trial logging and XDf data for a quick check for EML1

% synchronisation: absolute time on PC1 will be the master clock (to
% facilitate correspondance to the eyetracker events)

datapath = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EyeMindLink\Data';

sublist = [43];

for s = 1:length(sublist)
    pID = ['EML1_',sprintf('%03d',sublist(s))];
    % for trigger timing to nto get messed up, we need to disable jitter
    % removal on xdf load, otherwise it makes itt a regular stream.
    xdf_uncorr = load_xdf(fullfile(datapath, pID, [pID,'.xdf']),'HandleJitterRemoval' , false,'Verbose',true,'HandleClockSynchronization',false);
    xdf_corr = load_xdf(fullfile(datapath, pID, [pID,'.xdf']),'HandleJitterRemoval' ,false,'Verbose',true,'HandleClockSynchronization',true);
    
    % detect which stream is the trigger
    trig_ix = find(strcmp(cellfun(@(sas) sas.info.name, xdf_uncorr,'uni',false),{'eyeLink_trigger'}));
    xdf_trig_pc1 = xdf_corr{trig_ix}; 
    xdf_trig_pc2 = xdf_uncorr{trig_ix};  
    
    
    %% get eeg data start time in absolute time, taken from PC2 clock
    % note this value is found in the coresponding .vmrk file, todo write
    % code to pull this automatically!
    temp = fileread(fullfile(datapath, pID, 'EEG', [pID, '.vmrk']));
    eeg_start = regexp(temp, 'Mk1=New Segment,,1,1,0,(?<tst>[0-9]{20})','names').tst;
    eeg_start = [eeg_start(1:4) '-' eeg_start(5:6) '-' eeg_start(7:8) ' ' eeg_start(9:10) ':' eeg_start(11:12) ':' eeg_start(13:14) '.' eeg_start(15:end)];
    eeg_start_pc2abs = datetime(eeg_start, 'Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
    
    %% get the clock offset from the xdf
    % this is offset between boot times not between absolute times.
    % there are multiple measurements of the offset = we need to decide how to deal w that!
    % - Simplest is to compute mean offsest, use as a constant, then consider finer adjustment later w the residuals
    % - clock_offsets are the same in teh corrected and uncorrected files -
    % its just part of the XDF header
    all_offset_tpc1 = str2double({cell2mat(xdf_trig_pc1.info.clock_offsets.offset).time})';
    all_offset_val = str2double({cell2mat(xdf_trig_pc1.info.clock_offsets.offset).value})'; % note this is pc2 time minus pc1
    
    mean_offset = mean(all_offset_val); % in sec
    range_offset = range(all_offset_val); % in sec
    resid_offset = all_offset_val - mean_offset; % in sec, will use later to fine tune timings
    figure(99); clf
    scatter(all_offset_tpc1,all_offset_val,'.')
    xlabel('time (s since PC2 boot)')
    ylabel('offset (boot time PC2 - PC1)')
title([pID ' LSL clock offsets'], 'interpreter','none')

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
    
    %% COmpare log and LSL timestamps
    figure(22);clf
    plot(logtrig.datetime(1:length(xdf_trig_pc1.time_stamps)),xdf_trig_pc2.time_stamps,'o' )
    title([pID ' triallog vs LSL triggers'], 'interpreter','none')
    xlabel('log timestamp')
    ylabel('LSL timestamp')
    
    %% Estimate PC1 absolute boot time
    % subtract the LSL triggers uncorrected timestamps from the log
    % timestamps, should all be the same or v similar
    % note there are some weird delays introduced later in the session
    % where lsl timestamps lag behind the log, so only take events upt o
    % event_cutoff
    event_cutoff = 30;
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
    
    % % 1. earliest time
    %pc1_boot_abs = pc1_boot_ests(1);
    % 2. linear fit
    
    pc1_boot_fn = polyfit(xdf_trig_pc1.time_stamps(1:event_cutoff),seconds(pc1_boot_ests-min(pc1_boot_ests)),1);
    figure(303);hold on;
    x=linspace(0,max(xdf_trig_pc1.time_stamps),100);
    plot(x,seconds(x*pc1_boot_fn(1)+pc1_boot_fn(2))+min(pc1_boot_ests))
    pc1_boot_abs = pc1_boot_fn(2)+min(pc1_boot_ests);
    
%     pc1_boot_fn = polyfit(xdf_trig_pc1.time_stamps,datenum(pc1_boot_ests),1);
%     figure(303);hold on;
%     x=linspace(0,max(xdf_trig_pc1.time_stamps),5);
%     plot(x,datetime(x*pc1_boot_fn(1)+pc1_boot_fn(2),'ConvertFrom','datenum'))
%     pc1_boot_abs = datetime(pc1_boot_fn(2),'ConvertFrom','datenum','Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
    
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
    title([pID ' estimation of boot time PC2'],'interpreter','none')
    xlabel('timestamp (PC2)'); ylabel('timestamp (PC1)')
        %%    Load EEG using Fieldtrip
    % TODO: Is the ft function modifying timestamps to start at 0 or is
    % thsi in the recording?
    cfg=[];
    cfg.dataset = fullfile(datapath, pID, 'EEG', [pID, '.eeg']);
    eeg_bv = ft_preprocessing(cfg);
    %% Get EEG start time in other formats. Add this to the 0-based indexing of EEG timestamps to get it in that format
    eeg_start_pc2rel = seconds(eeg_start_pc2abs-pc2_boot_abs);
    eeg_start_pc1rel = offset_fn_tpc2(1)*eeg_start_pc2rel + offset_fn_tpc2(2);
    eeg_start_pc1abs = seconds(eeg_start_pc1rel)+pc1_boot_abs;
    eeg_start_pc1abs.Format=('yyyy-MM-dd HH:mm:ss.SSSSSS');
    
    %% Plot EEG and triggers from log file
    figure(s); clf;
    subplot(2,1,1)
    plot(seconds(eeg_bv.time{1,1})+eeg_start_pc1abs,eeg_bv.trial{1,1})
    hold on
    % plot vertical bars for the triggers
    plot([logtrig.datetime'; logtrig.datetime'],...
        repmat(ylim',1,length(logtrig.datetime)))
    xlabel('time')
    ylabel('raw EEG')
    subplot(2,1,2)
    plot(logtrig.datetime, logtrig.VAL, 's')
    xlabel('time')
    ylabel('trigger value')
    
    sgtitle([pID ' aligned to PC1 clock'], 'interpreter','none')
    
    %% store diagnistics
    timinfo(s).eeg_start_pc1abs = eeg_start_pc1abs;
    timinfo(s).eeg_start_pc2abs = eeg_start_pc2abs;
    timinfo(s).eeg_start_pc1rel = eeg_start_pc1rel;
    timinfo(s).eeg_start_pc2rel = eeg_start_pc2rel;
    timinfo(s).pc1_boot_abs = pc1_boot_abs;
    timinfo(s).pc2_boot_abs = pc2_boot_abs;
    timinfo(s).eeg_start_discrepancy = seconds(eeg_start_pc2abs - eeg_start_pc1abs);
    timinfo(s).pc1_boot_discrepancy = seconds(pc1_boot_ests(end) - pc1_boot_ests(1)); 
    timinfo(s).lsl_range_offset = range_offset; 

    
end

%%  unused crap

%%
%    %% get offset between pc1 boot (just running the windows comman,d this is shitty 1s resolution!) and EEG recording start (pc2)
%    % note this asssumes both clocks tell the same time!
%    % startup_time = '20200724083848';
%    % using online converter this is equiv to
%    pc1_boot_unix = 1595601528000;
%    unix_offset_ms = eeg_start_unix-pc1_boot_unix;
%
%% Plot synced data
%     figure(s); clf;
%     subplot(2,1,1)
%     plot(eeg_bv.time{1,1},eeg_bv.trial{1,1})
%     hold on
%     % plot vertical bars for the triggers
%     plot([xdf_trig.time_stamps - unix_offset_ms/1000;...
%         xdf_trig.time_stamps - unix_offset_ms/1000], ...
%         repmat(ylim',1,length(xdf_trig.time_stamps)))
%     xlabel('time (sec since EEG recording began)')
%     ylabel('raw EEG')
%     subplot(2,1,2)
% plot(xdf_trig.time_stamps - unix_offset_ms/1000, xdf_trig.time_series, 's')
%     xlabel('time (sec since EEG recording began)')
%     ylabel('trigger value')
%
% sgtitle([pID 'rough aligned by boot time (second resolution'], 'interpreter','none')