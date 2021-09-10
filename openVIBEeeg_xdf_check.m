% Load EEG and XDf data for a qick check for EML1 - EEG streamed with
% openVIBE

datapath = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EyeMindLink\Data';

sublist = {'020'};

for s = 1:length(sublist)
    pID = ['EML1_',sublist{s}];
    % for trigger timing to nto get messed up (loader will assume it's a regular stream) , we need to disable jitter
    % removal on xdf load.
    xdf_uncorr = load_xdf(fullfile(datapath, pID, [pID,'.xdf']),'HandleJitterRemoval' , false);
    xdf_corr = load_xdf(fullfile(datapath, pID, [pID,'.xdf']), 'HandleJitterRemoval' , false,'Verbose',true,'HandleClockSynchronization',true);
    
    trig_ix = find(strcmp(cellfun(@(sas) sas.info.name, xdf_uncorr,'uni',false),{'eyeLink_trigger'}));
    eeg_ix = find(contains(cellfun(@(sas) sas.info.name, xdf_corr,'uni',false),{'EEG'}));
    
    xdf_trig = xdf_uncorr{trig_ix};
    xdf_eeg = xdf_corr{eeg_ix};
    
%% get the absolute timestamps from the log file on PC1
    logtrig = readtable(fullfile(datapath,pID,[pID '_Trials.txt']) );
    % combine date and time to a datetime obj
    logtrig.Var1.Format = 'yyyy-MM-dd HH:mm:ss.SSSSSS';
    logtrig.datetime = logtrig.Var1 + logtrig.Var2;
    logtrig.Var1.Format = 'yyyy-MM-dd';

    %% Estimate PC1 absolute boot time
    % subtract the LSL triggers uncorrected timestamps from the log
    % timestamps, should all be the same or v similar
    %  first, find which events from the log correspond to LSL triggers (some are offline, meaning no recorded LSL triggers)
    log_trig_start = finddelay(logtrig.Var7, xdf_trig_pc1.time_series);
    pc1_boot_ests = logtrig.datetime(log_trig_start+(1:length(xdf_trig_pc1.time_stamps))) - ...
        seconds(xdf_trig_pc2.time_stamps)';
    
    % how do these estimates look over time?
    figure(303);clf
    scatter(xdf_trig_pc1.time_stamps,pc1_boot_ests)
    title('estimation of boot time PC1')
    xlabel('rel timestamp (PC1)'); ylabel('estimate boot time')
    
    %  Any variability in these estimates indicates the logging and LSL
    %  timestamp functions are not getting the same clocktime as one another.
    %
    %  which estimate should we use for PC1 boot time?
    %  -- I'd say the one derived from the earliest trigger event
    %  we could actually extrapolate this back in time to timestamp (PC1) ==0 if it looks like a
    %  simple function, or take an average over several early events
    
    pc1_boot_abs = pc1_boot_ests(1);
    %%
    figure(s); clf;
    subplot(2,1,1)
    % plot(eeg_bv.time{1,1},eeg_bv.trial{1,1})
    plot(xdf_eeg.time_stamps , xdf_eeg.time_series)
    hold on
    % plot vertical bars for the triggers
    plot([xdf_trig.time_stamps - str2double(xdf_eeg.info.created_at);...
        xdf_trig.time_stamps - str2double(xdf_eeg.info.created_at)], ...
        repmat(ylim',1,length(xdf_trig.time_stamps)))
    xlabel('time (sec since openVIBE stream began)')
    ylabel('raw EEG')
    subplot(2,1,2)
    plot(xdf_trig.time_stamps - str2double(xdf_eeg.info.created_at), xdf_trig.time_series, 's')
    xlabel('time (sec since openVIBE stream began)')
    ylabel('trigger value')
    
sgtitle(pID, 'interpreter','none')

%% get EEG header info 

%% process EEG in Fieldtrip

end