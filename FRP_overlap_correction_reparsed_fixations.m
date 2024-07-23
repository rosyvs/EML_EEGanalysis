% Use dataviewer reparsed and realigned fixations
% Extract condition-agnostic indivudual FRPs
% Save fixation-wise .set with info on whether blinks occurred in each epoch
% RVS 2024

clear all; close all force
init_unfold

% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
%%%%%%%%%
repro= 0; % re do analysis or just read in from file?
sublist = [34:181]; % TODO: replace with full sublist, short list used for dev
%%%%%%%%%

exclude_linenoise = [30 36 98 101 102 109 111 114 118 122 125 131 134 136 139]; % TODO: deal with this line noise
exclude_noise = [32 86]; %n other noise such as excessive jumps, blinks -  131? , movement
exclude_movement = [];
exclude_missingevents = [52 57 73 111 120 153]; % 57 73 because no eyetracking, others because EEG stopped recording early
exclude_noEEG = [1:18 23 77 88 138];
exclude_other = [22:24 26 27 31 39 40 78 160]; % TODO find reason for these - no .set why?
exclude = unique([exclude_linenoise exclude_noise exclude_movement exclude_missingevents exclude_noEEG exclude_other]); % Subj to exclude because no eeg or no trigger etc.
sublist = sublist(~ismember(sublist,exclude) );
%sublist = sublist( ismember(sublist,find(hasTriggerList.sdcard==1)));

%dir_raw = '/Volumes/Blue1TB/localEyeMindLink/Data';
dir_events = '~/Emotive Computing Dropbox/Rosy Southwell/EyeMindLink/Processed/events/';
dir_pre = fullfile('/Volumes/Blue1TB/EEG_processed') ; % prepro in MNE
dir_in = fullfile(dir_pre, 'preprocessed_set');
dir_fif = fullfile(dir_pre, 'preprocessed_fif');

unfdir = 'unfolded_FRP_reparsed';
mkdir( fullfile(dir_pre,unfdir))

beh_data = readtable('../../Data/EML1_page_level.csv');

%% load reparsed & IA-labelled fixations
fix_reparsed = readtable('/Volumes/Blue1TB/EyeMindLink/DataViewer/DataViewer_EML1/Output/FixationReport_14feb2023.txt');
cols = split(fix_reparsed.RECORDING_SESSION_LABEL,'-');
fix_reparsed.ParticipantID = cols(:,1);
clear cols
for s = 1:length(sublist)
    try
        tic;
        close all; clear EEGft logtrig fixations
        %  fileID = fopen('osc_win_log.txt', 'a');
    
        pID = ['EML1_',sprintf('%03d',sublist(s))];
        try
            EEG = pop_loadset(fullfile(dir_in, [pID '.set']));
        catch
            disp(['Failed to load preprocessed EEG file from ' fullfile(dir_in, [pID '.set'])])
            continue
        end
    
        nEEGsamples = EEG.pnts;
        events_orig = struct2table( EEG.event); % contains eyetracker events AND experiment events
        fixations = fix_reparsed(string(fix_reparsed.ParticipantID)==pID,:);
    
    
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
    
        % set up coarser task trialtypes
        logtrig.TrialType = repmat({'other'}, height(logtrig), 1);
        logtrig.TrialType(logtrig.VAL==7) = {'reading'};
        logtrig.TrialType(logtrig.VAL==20) = {'sham'};
        logtrig.TrialType(logtrig.VAL>=12 & logtrig.VAL <=15) = {'question'};
        logtrig.TrialType(logtrig.VAL>=2 & logtrig.VAL <=3) = {'recal'};
        logtrig.TrialType(logtrig.VAL>20 ) = {'localizer'};
        logtrig.TrialType(logtrig.VAL==25) = {'resting'};
    
    
        %% make new event structure forreading pages with label for identifier
        task_events = logtrig(:,{'EVENT','TrialType','eeg_use_sample','duration_sec'});
        task_events = renamevars(task_events,{'EVENT','TrialType','eeg_use_sample','duration_sec'},{'identifier','task','latency','duration'});
        % reduce options for task events to "reading" "sham" "qa"  "localizer" "resting" "other"
        task_events.duration = task_events.duration * EEG.srate; % into samples
        % task_events = task_events(strcmp(task_events.type,'reading'),:);
        task_events.urevent=NaN(height(task_events),1);
        task_events.type = repmat({'task_event'},height(task_events),1);
        task_events.page_fixation_ix = repmat(0, height(task_events),1); % column needed for consistency but doesn't mean anything for this event type
    
        %% make events for reparsed fixations
        % Note: eyetracker and EEG sample rate are both 1000Hz for EML. If this
        % werent the case you would need to do some conversions to match the
        % time units
        fixations = outerjoin(fixations,task_events, 'Type','Left', 'MergeKeys',true);
        % fixations = fixations(strcmp(fixations.type, 'reading'),:);
        fixations.duration = fixations.CURRENT_FIX_DURATION;
        fixations.latency = fixations.latency + fixations.CURRENT_FIX_START;
        fixations.type = repmat({'fix_R'}, height(fixations),1);
        % remove events outside of EEG recording
        fixations = fixations(fixations.latency>=0,:);
        fixations = fixations(fixations.latency+fixations.duration <= EEG.pnts,:);
        fixations.page_fixation_ix = fixations.CURRENT_FIX_INDEX; % needed for aligning to word level features

        %% select reading eye events from orig info
        eye_events = events_orig(contains(events_orig.type, 'fix_') | contains(events_orig.type, 'sac_')| contains(events_orig.type, 'blink_'),:);
        eye_events.task = repmat({'none'},height(eye_events),1);
        eye_events.identifier = repmat({'none'},height(eye_events),1);
        % eye_events.MW = repmat({'NA'},height(eye_events),1);
    
        % add task and MW label to eye events corresponding to reading only
        for e = 1:height(task_events) 
            sel = find((eye_events.latency >= task_events.latency(e)) & (eye_events.latency < task_events.latency(e)+task_events.duration(e) ));
            identifier = task_events.identifier(e);
            eye_events.task(sel) = repmat(task_events.task(e),length(sel),1);
            % eye_events.MW(sel) = repmat(task_events.MW(e),length(sel),1);
            eye_events.identifier(sel) = repmat(identifier,length(sel),1);
        end
        % remove NaN latency events
        eye_events = eye_events(~isnan(eye_events.latency),:);
        eye_events.page_fixation_ix = repmat(0, height(eye_events),1); % column needed for consistency but doesn't mean anything for this event type, only for labelled fixations

        %% compare reparsed and old fixations

        old_fixations = eye_events(strcmp(eye_events.type, 'fix_R'),:);
        old_fixations = old_fixations(:,{'identifier','type','latency','duration'});
        new_fixations = fixations(:,{'identifier','type','latency','duration'});
        new_fixations.type = repmat({'fix_R_reparsed'}, height(new_fixations),1);
        disp([num2str(height(old_fixations)) ' fixations before. ' num2str(height(new_fixations)) ' reparsed ' ])
        comp = sortrows([old_fixations; new_fixations], 'latency');
        new_fixations_final = max(fixations.latency)
        old_fixations_final = max(old_fixations.latency)
        old_fixations_keep = eye_events(strcmp(eye_events.type, 'fix_R')&eye_events.latency>new_fixations_final,:);
        
        %% make new eye_events using old blinks, saccs and new fixations (and old fixations for post-task where new fixation_report cuts off)
        % just use right eye blinks for now
        blinks = eye_events(~cellfun('isempty',strfind(eye_events.type, 'blink_R')),:);
        saccades = eye_events(~cellfun('isempty',strfind(eye_events.type, 'sac_R')),:);

        % columns to use 
        cols = {'type','latency','duration','task','identifier','page_fixation_ix'};
        events_new = [blinks(:,cols); fixations(:,cols); old_fixations_keep(:,cols); task_events(:,cols); saccades(:,cols)]; 
        events_new = sortrows(events_new,'latency');
        events_new.task= fillmissing(events_new.task,'constant','none');
        events_new.identifier= fillmissing(cellstr(events_new.identifier),'constant','none');
    
        %% put eye events into EEG event
        EEG.event = table2struct(events_new);
        EEG = eeg_checkset(EEG);
    
        %% Filter and resampke
        EEG = pop_eegfiltnew(EEG, 'locutoff',1);
        EEG = pop_eegfiltnew(EEG,'hicutoff',10,'filtorder' ,6600); % TODO read https://www.sciencedirect.com/science/article/abs/pii/S0165027014002866
        EEG = pop_resample( EEG, 40);
    
        %% design
        cfgDesign = [];
        cfgDesign.eventtypes = {'fix_R','blink_R','sac_R','task_event'}; % use right eye
        cfgDesign.formula = {'y ~ 1 + cat(task)','y ~ 1', 'y ~ 1', 'y ~ 1'}; 
    
        EEG = uf_designmat(EEG,cfgDesign);
    
        EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-0.2 1.5]);
        % %% artifacts
        % winrej = uf_continuousArtifactDetect(EEG,'amplitudeThreshold',10000);
        % EEG = uf_continuousArtifactExclude(EEG,struct('winrej',winrej));

        %% GLM fit - deconvolution
        EEG = uf_glmfit(EEG);
        % (strictly speaking optional, but recommended)
        uf_dc = uf_condense(EEG);
    
        %% GLM fit - no deconvolution
        EEG_epoch = uf_epoch(EEG,'timelimits',[-0.2 1.5]);
        EEG_epoch = uf_glmfit_nodc(EEG_epoch);
        uf_nodc = uf_condense(EEG_epoch);
    
        %% Plot single trials
        cfg=[];
        cfg.alignto = 'fix_R';
        cfg.split_by = 'task';
        cfg = [fieldnames(cfg),struct2cell(cfg)].';

        uf_erpimage( EEG ,'type','modelled','addResiduals',0,'channel',1,cfg{:}) ; 
        % addResiduals 
        %  cfg.addResiduals (integer 0-2):  default 0. 
%                                   If 1: Adds the overlap-including residuals.
%                                         That is: y_cont - X_dc*beta
%                                   If 2: Adds the corresponding residuals,
%                                         this will be identical to "cfg.method='raw'"
%                                         except when cfg.keep/cfg.remove
%                                         is used
        set(gcf, 'Position',[   616   707   806   310])
        saveas(gcf,fullfile(dir_pre,unfdir,[pID '_deconv_trialwise.png']) )

        uf_erpimage(EEG,'type','raw','channel',1, cfg{:})
        set(gcf, 'Position',[   616   707   806   310])
        saveas(gcf,fullfile(dir_pre,unfdir,[pID '_raw_trialwise.png']) )
    
        %% plot estimated ERP compare to no-deconvolution approach
        % without deconv
        g = uf_plotParam(uf_nodc,'channel',1,'deconv',0,'baseline',[-0.2 0],'add_intercept' ,1, 'plotSeparate','event');
    
        % set nice colors
        for gg = 1:length(g.results.geom_line_handle)
            g.results.geom_line_handle(gg).Color = [200 0 0]/255;
        end
    
        % with deconv
        g = uf_plotParam(uf_dc,'channel',1,'deconv',1,'baseline',[-0.2 0],'add_intercept' ,1,'gramm',g, 'plotSeparate','event','figure',0);
        % set nice colors
        for gg = 1:length(g.results.geom_line_handle)
            g.results.geom_line_handle(gg).Color = [0 200 0]/255;
        end
        set(gcf, 'Position',[   616   707   806   310])
    
        legend('no deconv','deconv')
        saveas(gcf,fullfile(dir_pre,unfdir,[pID '.png']) )
    
        %% save processed EEG
        pop_saveset(EEG, 'filename',pID, 'filepath',fullfile(dir_pre,unfdir), 'savemode','onefile')
    catch ME
        beep
        rethrow(ME)
        disp(['Fail for ' pID])
    end
end