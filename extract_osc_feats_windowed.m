% Extract spectral band features for each page
clear all; close all force
% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
%%%%%%%%%
repro= 0; % re do TF analysis or just read in from file?
sublist = [ 1:158];
%%%%%%%%%
exclude_linenoise = [30 36 98 101 102 109 111 114 118 122 125 131 134 136 139];
exclude_noise = [32 86]; %n other noise such as excessive jumps, blinks -  131? , movement
exclude_missingevents = [52 57 73 111 120 153]; % 57 73 because no eyetracking, others because EEG stopped recording early
exclude_noEEG = [1:18 23 77 88 138];
exclude_other = [20:26 27 31 39 40 78 160]; % TODO find reason for these - no .set why?
exclude = unique([exclude_linenoise exclude_noise exclude_missingevents exclude_noEEG exclude_other]); % Subj to exclude because no eeg or no trigger etc. 
sublist = sublist(~ismember(sublist,exclude) );
%sublist = sublist( ismember(sublist,find(hasTriggerList.sdcard==1)));

%dir_raw = '/Volumes/Blue1TB/localEyeMindLink/Data';
dir_events = '~/Dropbox (Emotive Computing)/EyeMindLink/Processed/events/';
dir_pre = fullfile('..','..','Data','EEG_processed') ;
dir_in = fullfile(dir_pre, 'preprocessed_set');
dir_fif = fullfile(dir_pre, 'preprocessed_fif');

oscdir = 'osc_windowed';

mkdir( dir_pre,oscdir)
texts = {'Bias','CausalClaims','Hypotheses','Validity','Variables'};
frange_labels = {'delta','theta','alpha','beta','gamma'};
franges = [2 4; 4 8; 8 13; 13 30; 30 70]; % note gamma cutoff of 70 - chosen because below 1/f knee
fratios = {2 3; 2 4; 4 {2 3}}; % numerator  denominator, by band index
fratio_labels = {'theta/alpha','theta/beta', 'engagementIndex'};

keep_channels = {'all'};
winlength = 2000; % ms window length
wininterval = 1000; % ms interval between intervals (if less than winlength, epochs overlap)
%%
T_all = [];
Tw_all = [];



if repro==1
    eeglab nogui % sets path defaults

    for s = 1:length(sublist)
        tic;
        close all; clear EEGft logtrig
        fileID = fopen('osc_win_log.txt', 'a');
        
        pID = ['EML1_',sprintf('%03d',sublist(s))];
        fprintf(fileID, pID);
        fprintf(fileID,'\n');
        %     %% Loading eyeCA cleaned data (Rosy is skeptical about that data)
        
        try
            EEG = pop_loadset(fullfile(dir_in, [pID '.set']));
        catch
            disp(['Failed to load preprocessed EEG file from ' fullfile(dir_in, [pID '.set'])])
            continue
        end
        EEGft = eeglab2fieldtrip(EEG,'raw','none');
        nEEGsamples = size(EEGft.time{1},2);
        events = struct2table( EEG.event); % contains eyetracker events AND experiment events
        % truncate events after end of EEG
        
        
        
        %% select channels
        eegchannels=ft_channelselection({keep_channels{:},'-x_dir','-y_dir','-z_dir'},EEGft.label);
        % remove xyz accelerometer
        cfg=[];
        cfg.channel= eegchannels;
        EEGft = ft_selectdata(cfg,EEGft);
        
        %% read metadata
        % Read info txt to determine whether EEG+triggers are from SD card
        % (default) or streamed (backup) or LSL (double backup!)
        fileinfo = readtxtfile(fullfile(dir_fif, [pID '-info.txt']));
        
        % read events.csv for triggers and descriptions
        logtrig = readtable(fullfile(dir_events,[pID '_events.csv']));
        % copy the correct EEGsample column for use depending on triginfo
        if contains(fileinfo, 'LA0','IgnoreCase',false)
            logtrig.eeg_use_sample = logtrig.eegSD_sample_est;
        else
            logtrig.eeg_use_sample = logtrig.eeg_sample_est;
        end
        
        %% Prepro for band features
        cfg=[];
        cfg.demean = 'yes';
        cfg.lpfilter = 'yes';
        cfg.lpfiltord = 1;
        cfg.lpfreq = 100; % antialiasing filter
        
        cfg.hpfilter = 'yes';
        cfg.hpfiltord = 2;
        cfg.hpfreq = 0.1;
        EEGft = ft_preprocessing(cfg, EEGft);
        
        %% 
        cfg=[];
        ft_databrowser(cfg,EEGft)
        
        %% make new cfg.trl with page starts and ends
        % cfg.trl must be a MATRIX with column specification:
        % 1. start sample
        % 2. end sample
        % 3. offset (usually 0)
        % 4. event code
        % additional columns as desired
        trl=[];
        trl(:,1) = logtrig.eeg_use_sample;
        trl(:,2) = trl(:,1)+logtrig.duration_sec*EEGft.fsample;
        trl(:,3) = 0;
        trl(:,4) = logtrig.VAL;
        for i=1:length(texts)
            trl(contains(logtrig.Text, texts{i}),5) = i; % represent text
        end
        trl(:,6) = logtrig.PageNum; % page
        trl(:,7) = logtrig.duration_sec;
        
        % remove NaN rows
        trl = trl(~isnan(trl(:,1)),:);
        
        % select reading only (i.e. val = 7 or 20) and resting state (25)
        trl = trl(ismember(trl(:,4), [7 20 25] ), :);
        
        % if any trials are past end of EEG recording, truncate
        if any(any(trl(:,[1 2])>nEEGsamples))
            disp('events past end of EEG recording, Truncating trl.')
            trl=trl( trl(:,1)<=nEEGsamples & trl(:,2)<=nEEGsamples,:);
        end
        
        % if any trials have NEGATIVE values for estimated EEG sample, truncate
        if any(any(trl(:,[1 2])<1))
            disp('events precede start of EEG recording, Truncating trl.')
            trl=trl( trl(:,1)>0 & trl(:,2)>0,:);
        end
        
        %% expand trl to window each trial to 2-second windows
        trl2 = [];
        for i = 1:length(trl)
            trl_expanded =[];
            trl_expanded(:,1) = trl(i,1):wininterval:trl(i,2)-winlength; % start
            trl_expanded(:,2) = trl_expanded(:,1)+winlength-1; % end
            trl_expanded(:,3) = zeros(size(trl_expanded(:,2))); % offset
            trl_expanded(:,4) = repmat(trl(i,4),size(trl_expanded(:,2)));
            trl_expanded(:,5) = repmat(trl(i,5),size(trl_expanded(:,2)));
            trl_expanded(:,6) = repmat(trl(i,6),size(trl_expanded(:,2)));
            trl_expanded(:,7) = repmat(winlength/1000,size(trl_expanded(:,2)));%duration_sec
            trl_expanded(:,8) = 1:size(trl_expanded,1);% window index
            trl_expanded(:,9) =repmat(i,size(trl_expanded(:,2)));% orig trial index
            
            trl2 = [trl2; trl_expanded];
        end
        
        
        %% absolute values artefact detection
        % remove huge jumps before z-based thresholkding for blinks
        % to prevent biasing of threshold by extreme values
        cfg=[];
        cfg.trl=trl2;
        cfg.artfctdef.threshold.range=2000;
        cfg.artfctdef.threshold.channel     = {'AFF5h','AFF6h'}; % frontal channels
        cfg.artfctdef.threshold.bpfilter   = 'no';
        
        disp(['detecting JUMP artifacts for ' pID])
        [cfg, artifact_jump] = ft_artifact_threshold(cfg, EEGft);
        jrej=zeros(length(trl2),1);
        for e=1:length(trl2)
            
            % remove jump artefact trials
            if any(artifact_jump(:,2) > trl2(e,1) & artifact_jump(:,2) < trl2(e,2))
                jrej(e)=1;
            end
        end
        disp([num2str(sum(jrej)) '/' num2str(length(trl2)) ' windows contain large-magnitude jumps.'])
        fprintf(fileID,[num2str(sum(jrej)) '/' num2str(length(trl2)) ' windows contain large-magnitude jumps.']);
        
        trl2 = trl2(find(~jrej),:);
        
        % blink artefact detection
        cfg=[];
        cfg.trl=trl2;
        cfg.artfctdef.zvalue.channel     = {'AFF5h','AFF6h'}; % frontal channels
        cfg.artfctdef.zvalue.cutoff      = 3;
        cfg.artfctdef.zvalue.trlpadding  = 0;
        cfg.artfctdef.zvalue.artpadding  = 0.2;
        cfg.artfctdef.zvalue.fltpadding  = 1;
        cfg.artfctdef.zvalue.bpfilter   = 'yes';
        cfg.artfctdef.zvalue.bpfilttype = 'but';
        cfg.artfctdef.zvalue.bpfreq     = [1 15];
        cfg.artfctdef.zvalue.bpfiltord  = 4;
        cfg.artfctdef.zvalue.hilbert    = 'yes';
        cfg.artfctdef.zvalue.interactive = 'no';
        disp(['detecting BLINK artifacts for ' pID])
        [cfg, artifact_blink] = ft_artifact_zvalue(cfg, EEGft);
        
        %% use ET events to reject intervals contaminated with blinks
        % let's say -50 to +400 ms rel to a blink constitutes the blink
        % window, and reject any overlapping windows
        % note we don't have ET for resting state so the automatic
        % detection by z threshold is all we have
        blink_intervals = events(strcmpi(events.type,'blink_either_eye'),:);
        blink_intervals.blinkstart = blink_intervals.latency - .05* EEG.srate;
        blink_intervals.blinkend = blink_intervals.latency + .4* EEG.srate;
        
        erej1 = zeros(length(trl2),1); erej2 = zeros(length(trl2),1);
        for e=1:length(trl2)
            
            % reject on basis of ET blinks
            if any(blink_intervals.blinkend > trl2(e,1) & blink_intervals.blinkstart < trl2(e,2))
                erej1(e)=1;
            end
            
            % reject on basis of Fieldtrips automatic blink detection from
            % frontal channels
            if any(artifact_blink(:,2) > trl2(e,1) & artifact_blink(:,2) < trl2(e,2))
                erej2(e)=1;
            end
        end
        erej=(logical(erej1) | logical(erej2));
        
        
        
        disp([num2str(sum(erej1)) ' windows overlap with ET blinks'])
        disp([num2str(sum(erej2))  ' windows overlap with EOG-detected blinks'])
        disp(['EOG detection recall = ' num2str(sum(logical(erej1) & logical(erej2)) / sum(erej1)) ])
        disp(['Removing '  num2str(sum(erej)) '/' num2str(length(erej)) ' epochs'])
        fprintf(fileID,'\n');
        
        fprintf(fileID,[num2str(sum(erej1)) ' windows overlap with ET blinks']);
        fprintf(fileID,'\n');
        fprintf(fileID,[num2str(sum(erej2))  ' windows overlap with EOG-detected blinks']);
        fprintf(fileID,'\n');
        fprintf(fileID,['EOG detection recall = ' num2str(sum(logical(erej1) & logical(erej2)) / sum(erej1)) ]);
        fprintf(fileID,'\n');
        fprintf(fileID,['Removing '  num2str(sum(erej)) '/' num2str(length(erej)) ' epochs']);
        fprintf(fileID,'\n');
        
        if sum(erej) / length(erej) >.5
            beep
            pause(0.25)
            beep
            disp("***WARNING*** blinks contaminate over half of windows. You may wish to exclude this participant***")
            fprintf(fileID,"***WARNING*** blinks contaminate over half of windows. You may wish to exclude this participant***")
            fprintf(fileID,'\n');
        end
        
        trl2 = trl2(find(1-erej),:);
        
        
        
        
        
        
        %% epoch into sham/reading/resting state windowed trials
        cfg=[];
        cfg.trl = trl2;
        EEGft = ft_redefinetrial(cfg, EEGft);
        
        % center the trials
        cfg=[];
        cfg.demean = 'yes';
        EEGft = ft_preprocessing(cfg, EEGft);
        
        
        %% resample - TODO does this resample the trl?
        cfg=[];
        cfg.detrend = 'yes';
        cfg.resamplefs = 200;
        EEGft = ft_resampledata(cfg, EEGft);
        
        %% Store metadata
        FREQ = [];
        FREQ.pID = pID;
        FREQ.trl=trl;
        FREQ.trl2=trl2;
        
        %% Freq analysis
        % https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_freq/#part-ii-spectral-analysis-on-eeg-resting-state-data
        % Welsh's method
        cfg             = [];
        cfg.output      = 'pow';
        cfg.channel     = 'all';
        cfg.method      = 'mtmfft';
        cfg.taper       = 'hanning';
        cfg.keeptrials  = 'yes';
        cfg.avgoverchan = 'yes';
        cfg.foi         = 2:.5:70;
        
        FREQ.PSDwin      =ft_selectdata(cfg,ft_freqanalysis(cfg, EEGft)); % this pads time to maximum trial duration, and gives NaN where this time is not present in a trial.
        % caution! Check these NaN values don't influence feature computaiton
        
        FREQ.freq = FREQ.PSDwin.freq;
        
        %% get PSD average over all windows per trialtype
        cfg=[];
        cfg.avgoverchan = 'yes';
        cfg.avgoverrpt = 'yes';
        cfg.nanmean = 'yes';
        
        cfg.trials = FREQ.PSDwin.trialinfo(:,1)==7;
        FREQ.PSDreading = ft_selectdata(cfg,FREQ.PSDwin);
        
        cfg.trials = FREQ.PSDwin.trialinfo(:,1)==20;
        FREQ.PSDsham = ft_selectdata(cfg,FREQ.PSDwin);
        
        cfg.trials = FREQ.PSDwin.trialinfo(:,1)==25;
        if(sum(cfg.trials) > 0)
            hasrest=1;
        else
            hasrest=0;
        end
        FREQ.PSDrest = ft_selectdata(cfg,FREQ.PSDwin);
        
        h=figure(1); clf
        plot(FREQ.PSDreading.freq, FREQ.PSDreading.powspctrm, 'b')
        hold on
        plot(FREQ.PSDsham.freq, FREQ.PSDsham.powspctrm, 'r')
        hold on
        plot(FREQ.PSDrest.freq, FREQ.PSDrest.powspctrm, 'k')
        legend({'reading','sham','rest'})
        title([pID ' PSD reading vs sham'],'Interpreter','none')
        export_fig(fullfile(dir_pre, oscdir, [pID, '_pagePSD.png']),'-png' )
        
        %% initialise table for page-level and window-level features
        
        T = table(trl(:,4), trl(:,5), trl(:,6),'VariableNames',{'VAL','texti','page'});
        Tw = table(trl2(:,4), trl2(:,5), trl2(:,6), trl2(:,8), trl2(:,9),'VariableNames',{'VAL','texti','page','win_ix','trl_ix'});
        
        %% %% Band-wise features %% %%
        %% put window-wise band power in Tw table
        for f= 1:length(franges) %frequency bands
            % select this frequency range and return all windows of this trial, avergae over freq only
            cfg=[];
            cfg.frequency = franges(f,:);
            cfg.avgoverfreq = 'yes';
            cfg.nanmean = 'yes';
            cfg.feedback = 'no';
            this_band_win = ft_selectdata(cfg,FREQ.PSDwin);
            pow =  this_band_win.powspctrm;
            pownorm = pow*length(this_band_win.freq)./squeeze(sum(FREQ.PSDwin.powspctrm,3)); % divide by total power this window (integral of PSD)
            
            %   pow_manual =  sqrt(nanmean(squeeze(FREQ.PSDwin.powspctrm(:,:,FREQ.PSDwin.freq >=franges(f,1) & FREQ.PSDwin.freq <franges(f,2))).^2,2));
            
            Twi = array2table([log(pow),log(pownorm)],'VariableNames',{'pow','pownorm'});
            desc = [frange_labels{f} '_'];
            Twi.Properties.VariableNames = strcat(desc, Twi.Properties.VariableNames) ;
            Tw = [Tw Twi];
        end
        
        
        %% Summarise over windows per original trial
        % full range PSD
        % use fieldtrip descriptives rather than matlab mean, variance
        cfg=[];
        cfg.variance = 'yes';
        cfg.keeptrials = 'no';
        clear power_total
        for i = 1:length(trl)
            cfg.trials = FREQ.PSDwin.trialinfo(:,6) == i;
            cfg.feedback = 'no';
            thistrial = ft_freqdescriptives(cfg,FREQ.PSDwin);
            FREQ.PSDtrial_mean(i,:) = thistrial.powspctrm;
            FREQ.PSDtrial_sem(i,:) =  thistrial.powspctrmsem;
            FREQ.PSDtrial_pow(i) = sum(FREQ.PSDtrial_mean(i,:));
            power_total(i) = log(FREQ.PSDtrial_pow(i));
        end
        T=[T table(power_total','VariableNames',{'power_total'})];
        % frequency band averages
        clear pow powsem pownorm
        for f= 1:length(franges) %frequency bands
            for i = 1:length(trl)
                % select this frequency range and return all windows of this trial, avergae over freq only
                cfg=[];  % use fieldtrip descriptives rather than matlab mean, variance
                cfg.variance = 'yes';
                cfg.keeptrials = 'no';
                cfg.frequency = franges(f,:);
                cfg.avgoverfreq = 'yes';
                cfg.trials = FREQ.PSDwin.trialinfo(:,6) == i;
                cfg.nanmean = 'yes';
                cfg.feedback = 'no';
                thistrial_band_win = ft_selectdata(cfg,FREQ.PSDwin);
                
                % mean/variance over all windows of this trial
                cfg=[];
                cfg.variance = 'yes';
                cfg.keeptrials = 'no';
                cfg.feedback = 'no';
                thistrial_band = ft_freqdescriptives(cfg,thistrial_band_win);
                pow(i) = thistrial_band.powspctrm;
                powsem(i) = thistrial_band.powspctrm;
                
                
                % Band power as ratio of total power over full range
                % note that nuemrator requires multiplying by numebr of frequency
                % bands because pow is average not total power, whereas
                % power_trialwise is total
                % average_band / average-full range would not work because each
                % computed using a different denominator
                pownorm(i) = thistrial_band.powspctrm*length(thistrial_band.freq)./FREQ.PSDtrial_pow(i); % divide by total power (integral of PSD)
                
                
            end
            % put in subject table
            Ti = table(log(pow') , log(powsem') , log(pownorm'),'VariableNames',{'pow','powsem','pownorm'});
            desc = [frange_labels{f} '_'];
            Ti.Properties.VariableNames = strcat(desc, Ti.Properties.VariableNames) ;
            T = [T Ti];
            
        end
        
        
        
        
        %% %% Band ratios %% %%
        % compute on non-normalsied band powers
        % These have already been log transformed so subtraction is equiv to dvosion before log transform
        
        %% Window-wise
        for r=1:length(fratios)
            numerator = table2array(Tw(:,[frange_labels{fratios{r,1}} '_pow']));
            if ~iscell(fratios{r,2}) % for combined bands
                denominator =table2array(Tw(:,[frange_labels{fratios{r,2}} '_pow']));
            else
                bands_to_combine = cell2mat(fratios{r,2});
                denominator = zeros(height(Tw),1);
                for i = bands_to_combine
                    denominator = denominator+table2array(Tw(:,[frange_labels{i} '_pow']));
                end
            end
            Twi = table(numerator-denominator, 'VariableNames',fratio_labels(r));
            Tw = [Tw Twi];
        end
        %% Trialwise
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
            Ti = table(numerator-denominator, 'VariableNames',fratio_labels(r));
            T = [T Ti];
        end
        
        %% %% peak alpha freq %% %%
        % get PSD in band for each freq - average over time AND channels
        alf = franges(contains(frange_labels,'alpha'),:);
        
        cfg=[];
        cfg.frequency = franges(contains(frange_labels,'alpha'),:);
        cfg.nanmean = 'yes';
        
        %% for windows
        alpha_i = ft_selectdata(cfg, FREQ.PSDwin);
        
        alpha_flat=zeros(length(alpha_i.freq), size(alpha_i.powspctrm,1) );
        % remove log-log linear trend on each trial
        aperiodic_alphai_fit = [];
        for t=1:size(alpha_i.powspctrm,1)
            aperiodic_alphai_fit(t,:) = polyfit(log(alpha_i.freq), squeeze(log(alpha_i.powspctrm(t,:,:))), 1);% normailise by 1/f
            %         alpha_flat = detrend(log(squeeze(alpha_i.powspctrm)'));% normailise by 1/f, detrend is column-wise
            alpha_flat(:,t) = (squeeze(log(alpha_i.powspctrm(t,:,:))) - aperiodic_alphai_fit(t,1)*log(alpha_i.freq)' -  aperiodic_alphai_fit(t,2));
            
        end
        [am,ai]=max(alpha_flat);% then find peak
        alpha_peak = alpha_i.freq(ai)';
        figure(4);
        plot(alpha_i.freq, alpha_flat)
        title('Alpha range per window, 1/f removed')
        Tw = [Tw table(alpha_peak)];
        
        %% for trials - get modal alpha freq\
        clear alpha_peak
        for i = 1:length(trl)
            trials = FREQ.PSDwin.trialinfo(:,6) == i;
            alpha_peak(i) = mode(Tw.alpha_peak(trials));
        end
        T = [T table(alpha_peak','VariableNames',{'alpha_peak'})];
        
        %% %% aperiodic decomposition %%
        %% window-wise
        if hasrest
            aperiodic_fit_rest = polyfit(log( FREQ.PSDrest.freq), log( FREQ.PSDrest.powspctrm), 1);
            aperiodic_exponent_rest = aperiodic_fit_rest(1);
            aperiodic_DC_rest = aperiodic_fit_rest(2);
            
            figure(5);clf
            plot(log( FREQ.PSDrest.freq), log( FREQ.PSDrest.freq)*aperiodic_exponent_rest+aperiodic_DC_rest)
            hold on
            plot(log( FREQ.PSDrest.freq), log( FREQ.PSDrest.powspctrm),'r')
            title('1/f fit on resting state')
        else
            aperiodic_exponent_rest = NaN; aperiodic_DC_rest = NaN;
        end
        
        % Fit trialwise
        aperiodic_fit=zeros(length(trl2) , 2);
        aperiodic_exponent = []; aperiodic_DC=[];
        for t=1:length(trl2)
            aperiodic_fit(t,:) = polyfit(log(FREQ.freq), squeeze(log(FREQ.PSDwin.powspctrm(t,:,:))),1);
            aperiodic_exponent(t,1) = aperiodic_fit(t,1);
            aperiodic_DC(t,1) = aperiodic_fit(t,2);
        end
        
        % Simialrity to resting state spectrum
        aperiodic_exponent_rel2rest = aperiodic_exponent - aperiodic_exponent_rest; % subtract rahter than divide as we are in loglog space
        aperiodic_DC_rel2rest = aperiodic_DC - aperiodic_DC_rest; % subtract rahter than divide as we are in loglog space
        Tw = [Tw table(aperiodic_exponent) table(aperiodic_DC), table(aperiodic_exponent_rel2rest),table(aperiodic_DC_rel2rest)];
        
        
        
        %% Fit trialwise
        aperiodic_fit=zeros(length(FREQ.PSDtrial_mean) , 2);
        aperiodic_exponent = []; aperiodic_DC=[];
        for t=1:length(trl)
            aperiodic_fit(t,:) = polyfit(log(FREQ.freq), log(FREQ.PSDtrial_mean(t,:)),1);
            aperiodic_exponent(t,1) = aperiodic_fit(t,1);
            aperiodic_DC(t,1) = aperiodic_fit(t,2);
        end
        
        % Simialrity to resting state spectrum
        aperiodic_exponent_rel2rest = aperiodic_exponent - aperiodic_exponent_rest; % subtract rahter than divide as we are in loglog space
        aperiodic_DC_rel2rest = aperiodic_DC - aperiodic_DC_rest; % subtract rahter than divide as we are in loglog space
        T = [T table(aperiodic_exponent) table(aperiodic_DC), table(aperiodic_exponent_rel2rest),table(aperiodic_DC_rel2rest)];
        
        %% %% SSA - spectral similarity analysis - see Zhou
        %  window wise
        SSA=[];
        if hasrest
            for t=1:length(trl2)
                SSA(t,:) = corr(FREQ.PSDrest.powspctrm', squeeze(FREQ.PSDwin.powspctrm(t,:,:)));
            end
        else
            SSA = NaN(length(trl2),1);
        end
        Tw = [Tw table(SSA) ];
        
        % whole trial
        SSA_wholetrial=[];
        if hasrest
            for t=1:length(trl)
                SSA_wholetrial(t,:) = corr(FREQ.PSDrest.powspctrm', FREQ.PSDtrial_mean(t,:)');
            end
        else
            SSA_wholetrial = NaN(length(trl),1);
        end
        T = [T table(SSA_wholetrial)];
        
        
        %% check correlation coefficient between features
        h= figure(3);
        colormap( getPyPlot_cMap('RdBu') )
        imagesc(corr(table2array(T(:,varfun(@isnumeric,T,'output','uniform'))),'Rows','pairwise'))                                 % Original 'XTick' Values
        xtlbl = T.Properties.VariableNames;                     % New 'XTickLabel' Vector
        set(gca, 'XTick',1:length(xtlbl), 'XTickLabel',xtlbl, 'XTickLabelRotation',90)
        set(gca, 'YTick',1:length(xtlbl), 'YTickLabel',xtlbl)
        set(gca,'TickLabelInterpreter','none')
        set(gcf,'Position',[1 1 800 800])
        export_fig(fullfile(dir_pre, oscdir, [pID, '_featcorrel.png']),'-png')
        
        h= figure(3);
        colormap( getPyPlot_cMap('RdBu') )
        imagesc(corr(table2array(Tw(:,varfun(@isnumeric,Tw,'output','uniform'))),'Rows','pairwise'))                                 % Original 'XTick' Values
        xtlbl = Tw.Properties.VariableNames;                     % New 'XTickLabel' Vector
        set(gca, 'XTick',1:length(xtlbl), 'XTickLabel',xtlbl, 'XTickLabelRotation',90)
        set(gca, 'YTick',1:length(xtlbl), 'YTickLabel',xtlbl)
        set(gca,'TickLabelInterpreter','none')
        set(gcf,'Position',[1 1 800 800])
        export_fig(fullfile(dir_pre, oscdir, [pID, '_winfeatcorrel.png']),'-png')
        
        %% table
        T.text(T.texti>0) = [texts(T.texti(T.texti>0))];
        %T_all.text(T_all.texti==0) = {'sham'};
        T.qType(T.VAL==7)={'reading'};
        T.qType(T.VAL==20)={'sham'};
        T.qType(T.VAL==25)={'restingState'};
        T.pID = repmat(pID, height(T),1);
        Tw.text(Tw.texti>0) = [texts(Tw.texti(Tw.texti>0))];
        %T_all.text(T_all.texti==0) = {'sham'};
        Tw.qType(Tw.VAL==7)={'reading'};
        Tw.qType(Tw.VAL==20)={'sham'};
        Tw.qType(Tw.VAL==25)={'restingState'};
        Tw.pID = repmat(pID, height(Tw),1);
        
        
        %% save TF decomposed data
        FREQ.pagefeat_table = T;
        FREQ.winfeat_table = Tw;
        FREQ.pID = pID;
        FREQ.franges = franges;
        FREQ.dir_in=dir_in;
        
        
        
        
        save(fullfile(dir_pre, oscdir,[pID '_FREQ.mat']),'FREQ');
        
        
        
        toc
        T_all = [T_all; T];
        Tw_all = [Tw_all; T];
        fprintf(fileID,'\n\n');
        
    end
fclose(fileID);
end

T_all=[]; Tw_all=[];
for s = 1:length(sublist)
    pID = ['EML1_',sprintf('%03d',sublist(s))];
    load(fullfile(dir_pre, oscdir,[pID '_FREQ.mat']),'FREQ');
    
    T_all = [T_all; FREQ.pagefeat_table];
    
    Tw_all = [Tw_all; FREQ.winfeat_table];
    
end
writetable(T_all,fullfile(dir_pre, oscdir, ['page_band_features_n' num2str(length(unique(cellstr(T_all.pID)))) '_upto_' pID '.csv']))
writetable(Tw_all,fullfile(dir_pre, oscdir, ['win_band_features_n' num2str(length(unique(cellstr(T_all.pID)))) '_upto_' pID '.csv']))
