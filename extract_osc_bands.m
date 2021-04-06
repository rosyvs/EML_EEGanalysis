% Extract spectral band features for each page
clear all; close all

% use only file w reliable trigger
triggerSources =readtable('triggerSources.csv');
sublist = find(triggerSources.sdcard+triggerSources.streamed + triggerSources.xdf); % subjects withat least one source of triggers
eeg_exclude = [20,21, 26, 73, 76, 77, 78, 99, 102]; % Subj to exclude because no eeg or no trigger or other issues
sublist = sublist(~ismember(sublist,eeg_exclude));
%sublist = find(hasTriggerList.sdcard==1);
%sublist = sublist(sublist~=73); % EEG but no eyetracking for these subjects

dir_raw = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EyeMindLink\Data\';
dir_info = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EML Rosy\Data\EEG_processed\';
dir_pre = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EML Rosy\Data\EEG_processed\';

mkdir( dir_pre, 'osc')
texts = {'Bias','CausalClaims','Hypotheses','Validity','Variables'};
frange_labels = {'delta','theta','alpha','beta','lowgamma','highgamma'};
franges = [2 4; 4 8; 8 13; 13 30; 30 60; 60 90];
fratios = {2 3; 2 4; 4 {2 3}}; % numerator  denominator, by band index
fratio_labels = {'theta/alpha','theta/beta', 'engagementIndex'}
cgroup_labels = {'frontocentral','occipitoparietal'};
cgroups = {{'CPz','FCz','AFF5h','AFF6h'},{'CCP5h','CCP6h','PPO9h','PPO10h'}};
%%
T_all = [];




for s = 1:length(sublist)
    tic;
    close all; clear EEGft logtrig
    pID = ['EML1_',sprintf('%03d',sublist(s))];
    
    %     %% Loading eyeCA cleaned data (Rosy is skeptical about that data)
    %eeglab nogui % sets path defaults
    
    %     dir_pre = 'C:\Users\roso8920\Dropbox (Emotive Computing)\EML Rosy\Data\EEG_processed\opticat_cleaned\';
    %
    %     EEGica = pop_loadset(fullfile(dir_pre, [pID '.set']));
    %     EEGft = eeglab2fieldtrip(EEGica,'raw','none');
    
    %% Loading basic preprocessed data from mne-python
    cfg = [];
    cfg.dataset = fullfile(dir_pre, [pID '_p.fif']);
    EEGft = ft_preprocessing(cfg);
    % note that MNE-py uses Volt but Fieldtrip generally uses uV.
    EEGft.trial{1} = EEGft.trial{1}*1000000;
    EEGft.elec.chanunit=repmat({'uV'},length(EEGft.elec.chanunit),1);
    
    % eeg channels
    eegchannels=ft_channelselection({'all','-x_dir','-y_dir','-z_dir'},EEGft.label);
    % remove xyz accelerometer
    cfg=[];
    cfg.channel= eegchannels;
    EEGft = ft_selectdata(cfg,EEGft);
    
    %% read logtrig
    % Read info txt to determine whether EEG+triggers are from SD card
    % (default) or streamed (backup)
    fileinfo = readtxtfile(fullfile(dir_info, [pID '-info.txt']));
    
    % read events.csv for triggers and descriptions
    logtrig = readtable(fullfile(dir_raw,pID,'EEG','events.csv'));
    % copy the correct EEGsample column for use depending on triginfo
    if contains(fileinfo, 'LA0','IgnoreCase',false)
        logtrig.eeg_use_sample = logtrig.eegSD_sample_est;
    else
        logtrig.eeg_use_sample = logtrig.eeg_sample_est;
    end
    
    %% make new cfg.trl with page starts and ends
    % cfg.trl must be a MATRIX with column specification:
    % 1. start sample
    % 2. end sample
    % 3. offset (usually 0)
    % 4. event code
    % additional columns as desired
    trl=[];
    trl(:,1) = logtrig.eeg_use_sample;
    % for trial duration, sometimes using responseTime.sec causes a small
    % overlap with the subsequent page. To prevent errors later (fieldtrip
    % doesn't like the overlapping trials) we take the start of the next
    % trial if the end time overlaps.
    logtrig.samp_to_next = diff([logtrig.eeg_use_sample; NaN]);
    durations=  nanmin(logtrig.responseTime_sec*EEGft.fsample , logtrig.samp_to_next-1);
    
    trl(:,2) = logtrig.eeg_use_sample + durations;
    
    trl(:,3) = 0;
    trl(:,4) = logtrig.VAL;
    for i=1:length(texts)
        trl(contains(logtrig.text, texts{i}),5) = i; % represent text
    end
    trl(:,6) = logtrig.page; % page
    
    % remove NaN rows
    trl = trl(~isnan(trl(:,1)),:);
    % select reading only (i.e. val = 7 or 20)
    trl = trl(ismember(trl(:,4), [7 20] ), :);
    
    % if any trials are past end of EEG recording, truncate
    if any(trl(:,[1 2])>size(EEGft.trial{1},2))
        disp('events past end of EEG recording, Truncating trl.')
        trl=trl( trl(:,1)<=size(EEGft.trial{1},2) & trl(:,2)<=size(EEGft.trial{1},2),:);
    end
    
    
    %% Prepro for band features
    cfg=[];
    cfg.demean = 'yes';
    cfg.lpfilter = 'yes';
    cfg.lpfiltord = 2;
    cfg.lpfreq = 100; % antialiasing filter
    
    cfg.hpfilter = 'yes';
    cfg.hpfiltord = 2;
    cfg.hpfreq = 0.1;
    cfg.reref  = 'yes';
    cfg.refchannel = 'all';
    cfg.refmethod = 'avg';
    EEGft = ft_preprocessing(cfg, EEGft);
    
    %% epoch into sham/reading trials
    cfg=[];
    cfg.trl = trl;
    EEGft = ft_redefinetrial(cfg, EEGft);
    
    % center the trials before detecting artefacts
    cfg=[];
    cfg.demean = 'yes';
    EEGft = ft_preprocessing(cfg, EEGft);
    
    
    %     %% reject segments of trials with artefacts
    %     cfg=[];
    %     %cfg.trl=trl;
    %     cfg.artfctdef.threshold.min   =-800;
    %     cfg.artfctdef.threshold.max   =800;
    %     cfg.artfctdef.threshold.bpfilter = 'no';
    %     cfg.continuous =  'no';
    %     [cfg,artefact] = ft_artifact_threshold(cfg,EEGft);
    %
    %
    %     cfg=[];
    %     cfg.artfctdef.reject = 'value'; % remove only the affected portion
    %     cfg.artfctdef.thresh.artifact = artefact; %
    %     cfg.artfctdef.value = 0;
    %     cfg.artfctdef.feedback        = 'no'
    %     EEGft = ft_rejectartifact(cfg,EEGft);
    % TODO remove known blinks
    
    
    %% resample
    cfg=[];
    cfg.detrend = 'yes';
    cfg.resamplefs = 200;
    EEGft = ft_resampledata(cfg, EEGft);
    
    %% TF analysis
    cfg            = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.keeptrials  = 'yes';
    cfg.pad = 'nextpow2';
    cfg.toi        = '50%'; % uses all timepoints, percentage specifies window overlap
    cfg.foi        = 2:0.1:90;
    cfg.t_ftimwin  = ones(size(cfg.foi)) * 2; % time window is the same for all freqa
    
    tfr      = ft_freqanalysis(cfg, EEGft);
    
    %% plot PSD overall for reading and sham over all channels
    cfg=[];
    cfg.avgoverchan = 'yes';
    cfg.avgovertime = 'yes';
    cfg.avgoverrpt = 'yes';
    cfg.nanmean = 'yes';
    
    PSDall = ft_selectdata(cfg,tfr); % average PSD for normalising each trial later
    cfg.trials = tfr.trialinfo(:,1)==7;
    PSDreading = ft_selectdata(cfg,tfr);
    cfg.trials = tfr.trialinfo(:,1)==20;
    PSDsham = ft_selectdata(cfg,tfr);
    
    h=figure(1); clf
    plot(PSDreading.freq, PSDreading.powspctrm, 'b')
    hold on
    plot(PSDsham.freq, PSDsham.powspctrm, 'r')
    legend({'reading','sham'})
    title([pID ' PSD reading vs sham'],'Interpreter','none')
    export_fig(fullfile(dir_pre, 'osc', [pID, '_pagePSD.png']),'-png' )
    
    %% average and get features
    h=figure(2); pc=0;
    T = table(tfr.trialinfo(:,1), tfr.trialinfo(:,2), tfr.trialinfo(:,3),'VariableNames',{'VAL','texti','page'});
    for f= 1:length(franges) %frequency bands
        pc=pc+1;
        
        % select this frequency range, avergae over freq only
        cfg=[];
        cfg.channel='all';
        cfg.frequency = franges(f,:);
        cfg.avgoverfreq = 'yes';
                    cfg.avgoverchan = 'yes';
        cfg.nanmean = 'yes';
        tfsel = ft_selectdata(cfg,tfr);
     
        subplot(length(franges),1, pc)
        plot(tfsel.time, squeeze(tfsel.powspctrm))
                title([frange_labels{f} '_power'],'Interpreter','none')

        cv = nanstd(squeeze(tfsel.powspctrm),0,2);
      sgtitle('Band power over time for each trial')
        % get PSD in band - average over time, but keep trials
        cfg=[];
        cfg.channel=eegchannels;
        cfg.frequency = franges(f,:);
        cfg.avgoverfreq = 'yes';
        cfg.avgoverchan = 'yes';
        cfg.avgovertime = 'yes';
        cfg.nanmean = 'yes';
        psdsel = ft_selectdata(cfg, tfr);

        % get PSD in band - average over time AND trials
        cfg=[];
        cfg.channel=eegchannels;
        cfg.frequency = franges(f,:);
        cfg.avgoverfreq = 'yes';
        cfg.avgoverchan = 'yes';
        cfg.avgovertime = 'yes';
        cfg.avgoverrpt = 'yes';
        cfg.nanmean = 'yes';
        avgpsdsel = ft_selectdata(cfg, tfr);
        
        pow = psdsel.powspctrm();
        cv = cv./pow;
        pownorm = pow./avgpsdsel.powspctrm; % divide by avergae over this subjects trials in this band
        
        % log transform power features
        pow=log(pow);
        pownorm=log(pownorm);
        
        % put in subject table
        Ti = table( pow, cv, pownorm);
        desc = [frange_labels{f} '_'];
            Ti.Properties.VariableNames = strcat(desc, Ti.Properties.VariableNames) ;
        T = [T Ti];
    end
    export_fig(fullfile(dir_pre, 'osc', [pID, '_bandpower.png']),'-png')
    
    %% Band ratios
    % compute on non-normalsied band powers
   % These have already been log transformed so subtraction is equiv to dvosion before log transform
    for r=1:length(fratios)
        numerator = table2array(T(:,[frange_labels{fratios{r,1}} '_pow']));
        if ~iscell(fratios{r,2}) % for combined bands
                            denominator =table2array(T(:,[frange_labels{fratios{r,2}} '_pow']));
        else
            bands_to_combine = cell2mat(fratios{r,2});
            denominator = zeros(height(T),1);
            for i = bands_to_combine
            denominator = denominator+table2array(T(:,[frange_labels{i} '_pow']));
            end
        end
        Tir = table(numerator-denominator, 'VariableNames',fratio_labels(r));
        T = [T Tir];
    end
    
    %% TODO peak alpha freq
     % get PSD in band for each freq - average over time AND channels
        cfg=[];
        cfg.channel = eegchannels;
        cfg.frequency = franges(contains(frange_labels,'alpha'),:);
        cfg.avgoverchan = 'yes';
        cfg.avgovertime = 'yes';
        cfg.nanmean = 'yes';
        alpha_i = ft_selectdata(cfg, tfr);
        [am,ai]=max(squeeze(alpha_i.powspctrm)');
    alpha_peak = alpha_i.freq(ai)';
    T = [T table(alpha_peak)];
    %% TODO aperiodic decomposition
    %% check correlation coefficient between features
    h= figure(3);
    imagesc(corrcoef(table2array(T)))                                 % Original 'XTick' Values
    xtlbl = T.Properties.VariableNames;                     % New 'XTickLabel' Vector
    set(gca, 'XTick',1:length(xtlbl), 'XTickLabel',xtlbl, 'XTickLabelRotation',90)
    set(gca, 'YTick',1:length(xtlbl), 'YTickLabel',xtlbl)
    set(gca,'TickLabelInterpreter','none')
    set(gcf,'Position',[1 1 800 800])
    export_fig(fullfile(dir_pre, 'osc', [pID, '_featcorrel.png']),'-png')
    
    toc
    T.pID = repmat(pID, height(T),1);
    T_all = [T_all; T];
end
T_all.text(T_all.texti>0) = [texts(T_all.texti(T_all.texti>0))];
T_all.text(T_all.texti==0) = {'sham'};

writetable(T_all,fullfile(dir_pre, 'osc', ['band_features_n' num2str(length(sublist)) '.csv']))

