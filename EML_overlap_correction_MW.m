% Extract spectral band features for each page
clear all; close all force
init_unfold

% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
%%%%%%%%%
repro= 0; % re do analysis or just read in from file?
sublist = [60 66 67 89 90];
%%%%%%%%%
exclude = [20:27 31 39 40 78 77 88 138 160]; % Subj to exclude because no eeg or no trigger etc.

exclude_linenoise = [30 36 98 101 102 109 111 114 118 122 125 131 134 136 139];
exclude_noise = [32 86]; %n other noise such as excessive jumps, blinks -  131? , movement
exclude_linenoise = [];
exclude_movement = [];

exclude_missingevents = [52 57 73 111 120 153]; % 57 73 because no eyetracking, others because EEG stopped recording early
exclude_noEEG = [1:18 23 77 88 138];
exclude_other = [22:24 26 27 31 39 40 78 160]; % TODO find reason for these - no .set why?
exclude = unique([exclude_linenoise exclude_noise exclude_missingevents exclude_noEEG exclude_other]); % Subj to exclude because no eeg or no trigger etc.
sublist = sublist(~ismember(sublist,exclude) );
%sublist = sublist( ismember(sublist,find(hasTriggerList.sdcard==1)));

%dir_raw = '/Volumes/Blue1TB/localEyeMindLink/Data';
dir_events = '~/Dropbox (Emotive Computing)/EyeMindLink/Processed/events/';
dir_pre = fullfile('..','..','Data','EEG_processed') ;
dir_in = fullfile(dir_pre, 'preprocessed_set');
dir_fif = fullfile(dir_pre, 'preprocessed_fif');

unfdir = 'unfolded_MW';

mkdir( dir_pre,unfdir)

beh_data = readtable('../../Data/EML1_page_level.csv');



for s = 1:length(sublist)
    tic;
    close all; clear EEGft logtrig
    %  fileID = fopen('osc_win_log.txt', 'a');

    pID = ['EML1_',sprintf('%03d',sublist(s))];
    try
        EEG = pop_loadset(fullfile(dir_in, [pID '.set']));
    catch
        disp(['Failed to load preprocessed EEG file from ' fullfile(dir_in, [pID '.set'])])
        continue
    end

    nEEGsamples = EEG.pnts;
    events = struct2table( EEG.event); % contains eyetracker events AND experiment events



    %% events
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

    % merge in behavioural data
    thisbeh = beh_data(string(beh_data.ParticipantID)==pID,:);
    thisbeh.TrialType = repmat({'reading'}, height(thisbeh), 1);
    logtrig = outerjoin(logtrig, thisbeh, 'Keys',{'TrialType', 'Text', 'PageNum'},'Type','left','MergeKeys',true);
    logtrig.MW(strcmp(logtrig.MW,'')) = {'NA'};
    % recode names to make NA the base case
    logtrig.MW(strcmp(logtrig.MW,'0')) = {'No'};
    logtrig.MW(strcmp(logtrig.MW,'1')) = {'Yes'};


    %% make new event structure column with label for reading vs sham
    task_events = logtrig(:,{'TrialType','eeg_use_sample','duration_sec','MW'});
    task_events = renamevars(task_events,{'TrialType','eeg_use_sample','duration_sec'},{'type','latency','duration'});
    task_events.duration = task_events.duration * EEG.srate; % into samples
    task_events = task_events(strcmp(task_events.type,'reading'),:);

    task_events.urevent=NaN(height(task_events),1);

    eye_events = events(contains(events.type, 'fix_') | contains(events.type, 'sac_')| contains(events.type, 'blink_'),:);
    eye_events.task = repmat({'none'},height(eye_events),1);
    eye_events.MW = repmat({'NA'},height(eye_events),1);

    for e = 1:height(task_events)
        sel = find((eye_events.latency >= task_events.latency(e)) & (eye_events.latency < task_events.latency(e)+task_events.duration(e) ));
        eye_events.task(sel) = repmat(task_events.type(e),length(sel),1);
        eye_events.MW(sel) = repmat(task_events.MW(e),length(sel),1);

    end


    % remove NaN latency events
    eye_events = eye_events(~isnan(eye_events.latency),:);

    EEG.event = table2struct(eye_events);
    EEG = eeg_checkset(EEG);

    %% Filter and resampke
    EEG = pop_eegfiltnew(EEG, 'locutoff',1);
    EEG = pop_eegfiltnew(EEG,'hicutoff',10,'filtorder' ,6600); % TODO read https://www.sciencedirect.com/science/article/abs/pii/S0165027014002866

    EEG = pop_resample( EEG, 40);

    % fieldtrip filter is way faster

    %                       EEGft = eeglab2fieldtrip(EEG,'raw','none');
    %
    %         cfg=[];
    %         cfg.demean = 'yes';
    %         cfg.lpfilter = 'yes';
    %         cfg.lpfiltord = 1;
    %         cfg.lpfreq = 100; % antialiasing filter
    %
    %         cfg.hpfilter = 'yes';
    %         cfg.hpfiltord = 2;
    %         cfg.hpfreq = 0.1;
    %         EEGft = ft_preprocessing(cfg, EEGft);
    %
    % clear header
    % header.Fs          =  EEG.srate;
    % header.nChans      =  EEG.nbchan;
    % header.nSamples    =  EEG.pnts;
    % header.nSamplesPre = -EEG.xmin*EEG.srate;
    % header.nTrials     =  EEG.trials;
    % header.label     = { EEG.chanlocs.labels }';
    %
    %   EEGfilt = fieldtrip2eeglab(header,EEGft);

    %% design
    cfgDesign = [];
    cfgDesign.eventtypes = {'fix_R'};
    % cfgDesign.formula = {'y ~ 1 + cat(task) + cat(MW)'};
    cfgDesign.formula = {'y ~ 1  + cat(MW)'};

    EEG = uf_designmat(EEG,cfgDesign);

    EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-0.2 1.5]);

    %% artifacts
    %     winrej = uf_continuousArtifactDetect(EEG,'amplitudeThreshold',1000);
    %     EEG = uf_continuousArtifactExclude(EEG,struct('winrej',winrej));
    %
    %% GLM fit - deconvolution
    EEG = uf_glmfit(EEG);
    % (strictly speaking optional, but recommended)
    uf_dc = uf_condense(EEG);

    %% GLM fit - no deconvolution
    EEG_epoch = uf_epoch(EEG,'timelimits',[-0.2 1.5]);

    EEG_epoch = uf_glmfit_nodc(EEG_epoch);
    uf_nodc = uf_condense(EEG_epoch);

    %% Plot single trials
    %     ax = uf_plotParam(uf_dc,'channel',2,'add_intercept' ,1);
    %
    %     ax2 = uf_plotParam(uf_nodc,'channel',1,'add_intercept' ,1);
    %
    uf_erpimage( EEG ,'type','modelled','addResiduals',1,'channel',1,'split_by','MW')
    %
    set(gcf, 'Position',[   616   707   806   310])
    saveas(gcf,fullfile(dir_pre,unfdir,[pID '_trialwise.png']) )
    % uf_erpimage(EEG,'type','raw','channel',1)

    %% plot compare to no-deconvolution approach
    % without deconv
    g = uf_plotParam(uf_nodc,'channel',1,'deconv',0,'baseline',[-0.2 0],'add_intercept' ,1);

    % set nice colors
    for gg = 1:length(g.results.geom_line_handle)
        g.results.geom_line_handle(gg).Color = [200 0 0]/255;
    end

    % with deconv
    g = uf_plotParam(uf_dc,'channel',1,'deconv',1,'baseline',[-0.2 0],'add_intercept' ,1,'gramm',g);
    % set nice colors
    for gg = 1:length(g.results.geom_line_handle)
        g.results.geom_line_handle(gg).Color = [0 200 0]/255;
    end

    set(gcf, 'Position',[   616   707   806   310])
    saveas(gcf,fullfile(dir_pre,unfdir,[pID '.png']) )

    legend('no deconv','deconv')
end