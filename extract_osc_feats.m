% Extract spectral band features for each page
clear all; close all
eeglab nogui % sets path defaults
ft_info off

% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
%%%%%%%%%
repro= 1; % re do TF analysis or just read in from file?
sublist = [ 121 160:181];
%%%%%%%%%
exclude_linenoise = [30 36 98 101 102 109 111 114 118 122 125 131 134 136 139];
exclude_noise = [32 86]; %n other noise such as excessive jumps, blinks -  131? , movement
exclude_missingevents = [52 57 73 111 120 153]; % 57 73 because no eyetracking, others because EEG stopped recording early
exclude_noEEG = [1:18 23 77 88 138 172];
exclude_other = [20:26 27 31 39 40 78 160 167 168 170 171 175 176 177 178 179]; % TODO find reason for these - no .set why? 170 onwards needs edf2asc.
exclude = unique([exclude_linenoise exclude_noise exclude_missingevents exclude_noEEG exclude_other]); % Subj to exclude because no eeg or no trigger etc.
sublist = sublist(~ismember(sublist,exclude) );
%sublist = sublist( ismember(sublist,find(hasTriggerList.sdcard==1)));

% dir_raw = '/Volumes/Blue1TB/localEyeMindLink/Data';
dir_events = '~/Dropbox (Emotive Computing)/EyeMindLink/Processed/events/';
dir_pre = fullfile('..','..','Data','EEG_processed') ;
dir_in = fullfile(dir_pre, 'preprocessed_set');
dir_fif = fullfile(dir_pre, 'preprocessed_fif');

oscdir = 'osc/'

mkdir( dir_pre,oscdir)
texts = {'Bias','CausalClaims','Hypotheses','Validity','Variables'};
frange_labels = {'delta','theta','alpha','beta','gamma'};
franges = [2 4; 4 8; 8 13; 13 30; 30 70]; % note gamma cutoff of 70 - chosen because below 1/f knee
fratios = {2 3; 2 4; 4 {2 3}}; % numerator  denominator, by band index
fratio_labels = {'theta/alpha','theta/beta', 'engagementIndex'};
cgroup_labels = {'frontocentral','occipitoparietal'};
cgroups = {{'CPz','FCz','AFF5h','AFF6h'},{'CCP5h','CCP6h','PPO9h','PPO10h'}};
%%
T_all = [];



if repro==1
    for s = 1:length(sublist)
        tic;
        close all; clear EEGft logtrig
        pID = ['EML1_',sprintf('%03d',sublist(s))];

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

        % eeg channels
        eegchannels=ft_channelselection({'all','-x_dir','-y_dir','-z_dir'},EEGft.label);
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

        %% use events to reject intervals with blinks



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

        %% epoch into sham/reading/resting state trials
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
        % https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_freq/#part-ii-spectral-analysis-on-eeg-resting-state-data
        % Welsh's method
        cfg             = [];
        cfg.output      = 'pow';
        cfg.channel     = 'all';
        cfg.method      = 'mtmconvol';
        cfg.taper       = 'hanning';
        cfg.keeptrials  = 'yes';
        cfg.pad         = 'nextpow2';
        cfg.toi         = '50%'; % uses all timepoints, percentage specifies window overlap
        cfg.foi         = 2:1:70;
        tfwinlength=2; % time window length in seconds
        cfg.t_ftimwin   = ones(size(cfg.foi)) * tfwinlength ; % time window is the same for all freqs

        TFR      = ft_freqanalysis(cfg, EEGft); % this pads time to maximum trial duration, and gives NaN where this time is not present in a trial.
        % caution! Check these NaN values don't influence feature computaiton

        %% PSD analysis for alpha - increased resolution in this range
        alf = franges(contains(frange_labels,'alpha'),:);
        cfg             = [];
        cfg.output      = 'pow';
        cfg.channel     = 'all';
        cfg.method      = 'mtmconvol';
        cfg.taper       = 'hanning';
        cfg.keeptrials  = 'yes';
        cfg.pad         = 'nextpow2';
        cfg.toi         = '50%'; % uses all timepoints, percentage specifies window overlap
        cfg.foi         = alf(1):0.1:alf(2);
        cfg.t_ftimwin   = ones(size(cfg.foi)) * tfwinlength; % time window is the same for all freqa

        TFRalpha      = ft_freqanalysis(cfg, EEGft); % this pads time to maximum trial duration, and gives NaN where this time is not present in a trial.
        % caution! Check these NaN values don't influence feature computaiton


        %% get PSD average over all channels and trials
        cfg=[];
        cfg.avgoverchan = 'yes';
        cfg.avgovertime = 'yes';
        cfg.avgoverrpt = 'yes';
        cfg.nanmean = 'yes';

        cfg.trials = TFR.trialinfo(:,1)==7;
        PSDreading = ft_selectdata(cfg,TFR);
        cfg.trials = TFR.trialinfo(:,1)==20;
        PSDsham = ft_selectdata(cfg,TFR);

        cfg.trials = TFR.trialinfo(:,1)==25;
        if(sum(cfg.trials) > 0)
            hasrest=1;
        else
            hasrest=0;
        end
        PSDrest = ft_selectdata(cfg,TFR);

        h=figure(1); clf
        plot(PSDreading.freq, PSDreading.powspctrm, 'b')
        hold on
        plot(PSDsham.freq, PSDsham.powspctrm, 'r')
        hold on
        plot(PSDrest.freq, PSDrest.powspctrm, 'k')
        legend({'reading','sham','rest'})
        title([pID ' PSD reading vs sham'],'Interpreter','none')
        export_fig(fullfile(dir_pre, oscdir, [pID, '_pagePSD.png']),'-png' )

        % Get trialwise PSD
        cfg=[];
        cfg.avgoverchan = 'yes';
        cfg.avgovertime = 'yes';
        cfg.nanmean = 'yes';
        cfg.avgoverrpt = 'no';
        PSDtrialwise = ft_selectdata(cfg,TFR);

        % get trialwise total power over all frequencies for normalized features
        power_trialwise = sum(PSDtrialwise.powspctrm,3);
        %% average and get features
        h=figure(2); pc=0;
        T = table(TFR.trialinfo(:,1), TFR.trialinfo(:,2), TFR.trialinfo(:,3),'VariableNames',{'VAL','texti','page'});
        for f= 1:length(franges) %frequency bands
            pc=pc+1;

            % select this frequency range, avergae over freq only
            cfg=[];
            cfg.channel='all';
            cfg.frequency = franges(f,:);
            cfg.avgoverfreq = 'yes';
            cfg.avgoverchan = 'yes';
            cfg.nanmean = 'yes';
            tfsel = ft_selectdata(cfg,TFR);

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
            psdsel = ft_selectdata(cfg, TFR);

            pow = psdsel.powspctrm();
            cv = cv./pow;

            % Band power as ratio of total power over full range
            % note that nuemrator requires multiplying by numebr of frequency
            % bands because pow is average not total power, whereas
            % power_trialwise is total
            % average_band / average-full range would not work because each
            % computed using a different denominator
            pownorm = pow*length(psdsel.freq)./power_trialwise; % divide by total power (integral of PSD)

            % !! Non-normalised power might be preferable, or normalising all
            % bands by the same denominator


            % log transform  features
            pow=log(pow);
            pownorm=log(pownorm);
            cv=log(cv);

            % put in subject table
            Ti = table( pow, cv, pownorm);
            desc = [frange_labels{f} '_'];
            Ti.Properties.VariableNames = strcat(desc, Ti.Properties.VariableNames) ;
            T = [T Ti];
        end
        export_fig(fullfile(dir_pre, oscdir, [pID, '_bandpower.png']),'-png')
        T=[T table(power_trialwise)];

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

        %% peak alpha freq
        % get PSD in band for each freq - average over time AND channels
        cfg=[];
        cfg.channel = eegchannels;
        cfg.frequency = franges(contains(frange_labels,'alpha'),:);
        cfg.avgoverchan = 'yes';
        cfg.avgovertime = 'yes';
        cfg.nanmean = 'yes';
        alpha_i = ft_selectdata(cfg, TFRalpha);

        alpha_flat=zeros(length(alpha_i.freq), size(alpha_i.powspctrm,1) );
        % remove log-log linear trend on each trial
        for t=1:size(alpha_i.powspctrm,1)
            aperiodic_alphai_fit(t,:) = polyfit(log(alpha_i.freq), squeeze(log(alpha_i.powspctrm(t,:,:))), 1);% normailise by 1/f
            %         alpha_flat = detrend(log(squeeze(alpha_i.powspctrm)'));% normailise by 1/f, detrend is column-wise
            alpha_flat(:,t) = (squeeze(log(alpha_i.powspctrm(t,:,:))) - aperiodic_alphai_fit(t,1)*log(alpha_i.freq)' -  aperiodic_alphai_fit(t,2));

        end
        [am,ai]=max(alpha_flat);% then find peak
        alpha_peak = alpha_i.freq(ai)';
        figure(4);
        plot(alpha_i.freq, alpha_flat)
        title('Alpha range, 1/f removed')
        T = [T table(alpha_peak)];

        %% aperiodic decomposition
        %  TODO run  on resting state?
        % this fit is on all reading/sham pages
        if hasrest
            aperiodic_fit_rest = polyfit(log(PSDrest.freq), log(PSDrest.powspctrm), 1);
            aperiodic_exponent_rest = aperiodic_fit_rest(1);
            aperiodic_DC_rest = aperiodic_fit_rest(2);


            figure(5);clf
            plot(log(PSDrest.freq), log(PSDrest.freq)*aperiodic_exponent_rest+aperiodic_DC_rest)
            hold on
            plot(log(PSDrest.freq), log(PSDrest.powspctrm),'r')
            title('1/f fit on resting state')
        else
            aperiodic_exponent_rest = NaN; aperiodic_DC_rest = NaN;
        end

        % this fit is trialwise
        aperiodic_fit=zeros(size(PSDtrialwise.powspctrm,1) , 2);
        aperiodic_exponent = []; aperiodic_DC=[];
        for t=1:size(PSDtrialwise.powspctrm,1)
            aperiodic_fit(t,:) = polyfit(log(PSDtrialwise.freq), squeeze(log(PSDtrialwise.powspctrm(t,:,:))), 1);
            aperiodic_exponent(t,1) = aperiodic_fit(t,1);
            aperiodic_DC(t,1) = aperiodic_fit(t,2);
        end
        T = [T table(aperiodic_exponent) table(aperiodic_DC)];

        %% Simialrity to resting state spectrum
        aperiodic_exponent_rel2rest = aperiodic_exponent - aperiodic_exponent_rest; % subtract rahter than divide as we are in loglog space
        aperiodic_DC_rel2rest = aperiodic_DC - aperiodic_DC_rest; % subtract rahter than divide as we are in loglog space
        %% SSA - spectral similarity analysis - see Zhou

        if hasrest
            for t=1:length(trl)
                SSA_wholetrial(t,:) = corr(PSDrest.powspctrm', squeeze(PSDtrialwise.powspctrm(t,1,:)));
            end
        else
            SSA_wholetrial = NaN(length(trl),1);
        end
        T = [T table(SSA_wholetrial)  table(aperiodic_exponent_rel2rest) table(aperiodic_DC_rel2rest)];
        %% check correlation coefficient between features
        h= figure(3);
        colormap( brewermap([],'RdBu') )
        imagesc(corr(table2array(T(:,varfun(@isnumeric,T,'output','uniform'))),'Rows','pairwise'))                                 % Original 'XTick' Values
        xtlbl = T.Properties.VariableNames;                     % New 'XTickLabel' Vector
        set(gca, 'XTick',1:length(xtlbl), 'XTickLabel',xtlbl, 'XTickLabelRotation',90)
        set(gca, 'YTick',1:length(xtlbl), 'YTickLabel',xtlbl)
        set(gca,'TickLabelInterpreter','none')
        set(gcf,'Position',[1 1 800 800])
        export_fig(fullfile(dir_pre, oscdir, [pID, '_featcorrel.png']),'-png')

        %% table
        T.text(T.texti>0) = [texts(T.texti(T.texti>0))];
        %T_all.text(T_all.texti==0) = {'sham'};
        T.qType(T.VAL==7)={'reading'};
        T.qType(T.VAL==20)={'sham'};
        T.qType(T.VAL==25)={'restingState'};
        T.pID = repmat(pID, height(T),1);


        %% save TF decomposed data
        TF_struct.feat_table = T;
        TF_struct.PSDtrialwise = PSDtrialwise;
        TF_struct.PSDreading = PSDreading;
        TF_struct.PSDsham = PSDsham;
        TF_struct.PSDrest = PSDsham;
        TF_struct.TFR = TFR;
        TF_struct.pID = pID;
        TF_struct.trl=trl;
        TF_struct.franges = franges;
        TF_struct.dir_in=dir_in;




        save(fullfile(dir_pre, oscdir,[pID '_TF_pages.mat']),'TF_struct');



        toc
        T_all = [T_all; T];

        clear SSA_wholetrial
    end
end

T_all=[];
for s = 1:length(sublist)
    pID = ['EML1_',sprintf('%03d',sublist(s))];
    load(fullfile(dir_pre, oscdir,[pID '_TF_pages.mat']),'TF_struct');

    T_all = [T_all; TF_struct.feat_table];


end
writetable(T_all,fullfile(dir_pre, oscdir, ['band_features_n' num2str(length(unique(cellstr(T_all.pID)))) '_upto_' pID '.csv']))
