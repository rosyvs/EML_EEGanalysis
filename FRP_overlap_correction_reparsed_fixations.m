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
sublist = [ 128:158]; % TODO: replace with full sublist, short list used for dev
%%%%%%%%%

exclude_linenoise = [30 36 98 101 102 109 111 114 118 122 125 131 134 136 139]; % TODO: deal with this line noise
exclude_noise = [32 86]; %n other noise such as excessive jumps, blinks -  131? , movement
exclude_movement = [];
exclude_missingevents = [52 57 73 111 120 153]; % 57 73 because no eyetracking, others because EEG stopped recording early
exclude_noEEG = [1:18 23 77 88 138 79 87 92 127]; % 79 onwards missing from reparsed fixations
exclude_other = [22:24 26 27 31 39 40 78 160]; % TODO find reason for these - no .set why?
exclude_glmfitfail = [19 35 68 147 149 ] ; % 54 onwards are only failing with sac splines
exclude = unique([exclude_linenoise exclude_noise exclude_movement exclude_missingevents exclude_noEEG exclude_glmfitfail exclude_other]); % Subj to exclude because no eeg or no trigger etc.
sublist = sublist(~ismember(sublist,exclude) );
%sublist = sublist( ismember(sublist,find(hasTriggerList.sdcard==1)));

%dir_raw = '/Volumes/Blue1TB/localEyeMindLink/Data';
dir_events = '~/Emotive Computing Dropbox/Rosy Southwell/EyeMindLink/Processed/events/';
dir_pre = fullfile('/Volumes/Blue1TB/EEG_processed') ; % prepro in MNE
dir_in = fullfile(dir_pre, 'preprocessed_set');
dir_fif = fullfile(dir_pre, 'preprocessed_fif');

unfdir = 'unfolded_FRP_reparsed_v6';
mkdir( fullfile(dir_pre,unfdir))

beh_data = readtable('../../Data/EML1_page_level.csv');
%% load reparsed & IA-labelled fixations
fix_reparsed =readtable('/Users/roso8920/Emotive Computing Dropbox/Rosy Southwell/EML Rosy/EMLLM/info/FixationReport+InboundSaccades.csv',"TextType","char");

cols = split(fix_reparsed.RECORDING_SESSION_LABEL,'-');
fix_reparsed.ParticipantID = cols(:,1);
clear cols
v=figure('Visible','off');
v.WindowState = 'minimized' ;
for s = 1:length(sublist)
        tic;
        clear EEGft logtrig fixations
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

        %% make events for reparsed fixations
        % Note: eyetracker and EEG sample rate are both 1000Hz for EML. If this
        % werent the case you would need to do some conversions to match the
        % time units
        fixations = outerjoin(fixations,task_events, 'Type','Left', 'MergeKeys',true);
        task_events.page_fixation_ix = repmat(0, height(task_events),1); % column needed for consistency but doesn't mean anything for this event type
        task_events.IA_ID = repmat({'.'}, height(task_events),1); % column needed for consistency but doesn't mean anything for this event type
        task_events.INBOUND_SAC_AMPLITUDE=NaN(height(task_events),1); %
        % fixations = fixations(strcmp(fixations.type, 'reading'),:);
        fixations.duration = fixations.CURRENT_FIX_DURATION;
        fixations.latency = fixations.latency + fixations.CURRENT_FIX_START;
        fixations.type = repmat({'fix_R'}, height(fixations),1);
        % remove events outside of EEG recording
        fixations = fixations(fixations.latency>=0,:);
        fixations = fixations(fixations.latency+fixations.duration <= EEG.pnts,:);
        fixations.page_fixation_ix = fixations.CURRENT_FIX_INDEX; % needed for aligning to word level features
        fixations.IA_ID = fixations.CURRENT_FIX_INTEREST_AREA_ID; % needed for aligning to word level features

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
        eye_events.IA_ID = repmat({'.'}, height(eye_events),1); % column needed for consistency but doesn't mean anything for this event type
        eye_events.INBOUND_SAC_AMPLITUDE=NaN(height(eye_events),1); %
        
        %% compare reparsed and old fixations
        old_fixations = eye_events(strcmp(eye_events.type, 'fix_R'),:);
        old_fixations = old_fixations(:,{'identifier','type','latency','duration'});
        new_fixations = fixations(:,{'identifier','type','latency','duration'});
        new_fixations.type = repmat({'fix_R_reparsed'}, height(new_fixations),1);
        disp([num2str(height(old_fixations)) ' fixations before. ' num2str(height(new_fixations)) ' reparsed ' ])
        comp = sortrows([old_fixations; new_fixations], 'latency');
        new_fixations_final = max(fixations.latency);
        old_fixations_final = max(old_fixations.latency);
        old_fixations_keep = eye_events(strcmp(eye_events.type, 'fix_R')&eye_events.latency>new_fixations_final,:);
        
        %% make new eye_events using old blinks, saccs and new fixations (and old fixations for post-task where new fixation_report cuts off)
        % just use right eye blinks for now
        blinks = eye_events(~cellfun('isempty',strfind(eye_events.type, 'blink_R')),:);
        saccades = eye_events(~cellfun('isempty',strfind(eye_events.type, 'sac_R')),:);

        % columns to use 
        cols = {'type','latency','duration','task','identifier','page_fixation_ix','IA_ID','INBOUND_SAC_AMPLITUDE'};
        events_new = [blinks(:,cols); fixations(:,cols); old_fixations_keep(:,cols); task_events(:,cols); saccades(:,cols)]; 
        events_new = sortrows(events_new,'latency');
        events_new.task= fillmissing(events_new.task,'constant','none');
        events_new.identifier= fillmissing(cellstr(events_new.identifier),'constant','none');
        clear fixations events

        %% put eye events into EEG event
        EEG.event = table2struct(events_new);
        EEG = eeg_checkset(EEG);
    
        %% trim data to just reading task (remove noisy setup)
        read_events = events_new(strcmp(events_new.task ,'reading'),:);
        margin = 60000; % milliseconds
        crop_at = [min(read_events.latency)-margin max(read_events.latency)+margin];
        EEG = pop_select(EEG, 'time', crop_at/1000); % contraty to docs it is in seconds

        %% impute missing INBOUND_SACCADE_AMP with median
        ev= struct2table(EEG.event);
        ev.INBOUND_SAC_AMPLITUDE{1}= NaN; % fill in for empty cells at "boundary" event
        ev.INBOUND_SAC_AMPLITUDE{end}= NaN; % fill in for empty cells at "boundary" event
        ev.INBOUND_SAC_AMPLITUDE = cell2mat(ev.INBOUND_SAC_AMPLITUDE); 

        tasks = unique({EEG.event.task});
        for tk = 1:length(tasks)
            task=tasks{tk};
            if (isempty(task) | strcmp(task,'boundary'))
                continue
            end
            md = median(ev.INBOUND_SAC_AMPLITUDE(strcmp(task,ev.task),:),'omitmissing');
          
            disp(['median sac amp for ' task ':' num2str(md)])
            ev.INBOUND_SAC_AMPLITUDE(strcmp(task,ev.task)) = fillmissing(ev.INBOUND_SAC_AMPLITUDE(strcmp(task,ev.task)),'constant',md);
        end

        EEG.event = table2struct(ev);
        %% Reref to average (might be better fro N400 as this should be maximal at central site but ref is at Cz)
        EEG = pop_reref(EEG,[],'interpchan',[]);

        %% Filter and resample (as per OG paper)
        EEG = pop_eegfiltnew(EEG, 'locutoff',.1,'usefftfilt',1); % as per unfold paper,
        % thouogh detrending should have fixed slow drifts & this is very
        % slow
        EEG = pop_eegfiltnew(EEG,'hicutoff',40,'usefftfilt',1); % TODO read https://www.sciencedirect.com/science/article/abs/pii/S0165027014002866
        sr_orig = EEG.srate;
        EEG = pop_resample( EEG, 100); % duration and latency fields of .event are now resampled too - untis are samples in the new srate
    
        %% design
        cfgDesign = [];
        cfgDesign.eventtypes = {'fix_R','task_event'}; % use right eye
        cfgDesign.formula = {'y ~ 1 + cat(task) + spl(INBOUND_SAC_AMPLITUDE,5)', 'y ~ 1'}; 
        % cfgDesign.eventtypes = {'fix_R','blink_R','sac_R','task_event'}; % use right eye
        % cfgDesign.formula = {'y ~ 1 + cat(task) + spl(CURRENT_SACCADE_AMPLITUDE,5)','y ~ 1+ spl(CURRENT_SACCADE_AMPLITUDE,5)', 'y ~ 1', 'y ~ 1'}; 
        % check what event types are in the data and omit any that dont exist otherwise it won;t make a design matrix
        valid_events = unique({EEG.event.type});
        keep_regressors = ones(length(cfgDesign.eventtypes),1);
        for j = 1:length(cfgDesign.eventtypes)
            if ~ismember(cfgDesign.eventtypes{j},valid_events)
                disp(['no isntances of ' cfgDesign.eventtypes{j}])
                keep_regressors(j)=0;
            end
        end
        cfgDesign.eventtypes = cfgDesign.eventtypes(find(keep_regressors));
        cfgDesign.formula = cfgDesign.formula(find(keep_regressors));
        EEG = uf_designmat(EEG,cfgDesign);
        if sum(sum(isnan(EEG.unfold.X))) >0
            disp(['Found ' mat2str(sum(isnan(EEG.unfold.X))) ' NaN values in columns of X. Imputing...']);
            EEG = uf_imputeMissing(EEG);
        end
        EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-0.3 0.8]);

        %% artifacts
        winrej = uf_continuousArtifactDetect(EEG,'amplitudeThreshold',[-1500 1500]);
        EEG = uf_continuousArtifactExclude(EEG,struct('winrej',winrej));

        %% GLM fit - deconvolution
        EEG = uf_glmfit(EEG);
        uf_dc = uf_condense(EEG); % (strictly speaking optional, but recommended)
    
        %% GLM fit - no deconvolution
        EEG_epoch = uf_epoch(EEG,'timelimits',[-0.3 0.8]);
        EEG_epoch = uf_glmfit_nodc(EEG_epoch);
        uf_nodc = uf_condense(EEG_epoch);

        %% Extract overlap-corrected single trial ERPs
        cfg=[];
        cfg.alignto = 'fix_R';
        % cfg.baseline = [-100, 0];
        cfg.winrej = winrej;
        cfg = [fieldnames(cfg),struct2cell(cfg)].';
        dc_FRP = uf_getERP(EEG ,'type','modelled','addResiduals',1,'channel',1:length(EEG.chanlocs),cfg{:});
        raw_FRP = uf_getERP(EEG ,'type','raw','addResiduals',0,'channel',1:length(EEG.chanlocs),cfg{:});
        dc_FRP_noresid = uf_getERP(EEG ,'type','modelled','addResiduals',0,'channel',1:length(EEG.chanlocs),cfg{:});

        %% plot FRPs and save single trials
        set(0,'CurrentFigure',v);clf
        colororder("glow12")
        channels = 1:length(dc_FRP.chanlocs);
        keep_evt = struct2table(dc_FRP.urevent);
        ax1=subplot(2,1,1);
        ax2=subplot(2,1,2);
        for c = channels
            dcFRPi = squeeze(dc_FRP.data(c,:,:));
            rawFRPi = squeeze(raw_FRP.data(c,:,:));

            % rawFRPi=squeeze(EEG_epoch.data(c,:,keep_evt.urevent));
            % checn n nunique values
            nuniqueEEG = length(unique(dcFRPi));
            % plot to check
            subplot(2,1,1)
            title('raw FRP')
            plot(dc_FRP.times,mean(rawFRPi,2),'LineWidth',2);
            hold on
            subplot(2,1,2)
            title('deconvolved FRP')
            plot(dc_FRP.times,mean(dcFRPi,2),'LineWidth',2);
            hold on
            chanlabel = dc_FRP.chanlocs(c).labels;
            writematrix(dcFRPi,fullfile(dir_pre,unfdir, [pID '_unfoldedFRP-reading_' chanlabel '.csv']));
            writematrix(rawFRPi,fullfile(dir_pre,unfdir, [pID '_rawFRP-reading_' chanlabel '.csv']));

        end
        sgtitle(pID)
        legend({dc_FRP.chanlocs.labels})
        saveas(gcf,fullfile(dir_pre,unfdir,[pID '_butterfly.png']) )

        %% N400 magnitude and latency
        % [mag, lat] = getERPmagnitude(dc_FRP, [300, 500], keep_ix
        win=[300 500];
        sel = find(dc_FRP.times<= win(2) & dc_FRP.times>=win(1));
        % what about winrej trials? are they a;ready removed? 
        for c = channels
            y=dc_FRP.data;
            y = y(c,sel,:);
            [m, ix] = max(abs(y),[],2);
            mag = zeros(size(ix(:))); % max val for each trial. I'm being dumb but cant work out how to vectorise this rn
            for z = 1:length(ix(:))% loop over trials
                mag(z) = y(1,ix(z),z);
            end
            lat = squeeze(dc_FRP.times(sel(squeeze(ix))));
            chanlabel = dc_FRP.chanlocs(c).labels;
            keep_evt.(['n400_magnitude_' chanlabel]) = mag;
            keep_evt.(['n400_latency_' chanlabel]) = lat';
        end
        % convert durations to ms from samples
        keep_evt.latency = keep_evt.latency*sr_orig/dc_FRP.srate;
        keep_evt.duration = keep_evt.duration*sr_orig/dc_FRP.srate;
        writetable(keep_evt, fullfile(dir_pre,unfdir,[pID '_reading_N400_stats.csv']));        %% Plot single trials
        
        %% plot single trials ERP image
        cfg=[];
        cfg.alignto = 'fix_R';
        cfg.winrej = winrej;
        cfg.sort_alignto = 'fix_R';
        cfg.sort_time = [eps 0.8];
        % % cfg.keep={{'task','reading'}}; % which event(s) predictor(s) pair(s) to keep   {'eventA',{'stimB'}} (or a cell array of such cell arrays)
        % cfg.split_by = 'task';
        cfg.figure=0;
        cfg = [fieldnames(cfg),struct2cell(cfg)].';

        set(0,'CurrentFigure',v); clf
        uf_erpimage(EEG,'type','raw','channel',1, cfg{:})
        % set(v, 'Position',[   616   707   806   310])
        saveas(gcf,fullfile(dir_pre,unfdir,[pID '_raw_trialwise.png']) )

        set(0,'CurrentFigure',v); clf
        xx=uf_erpimage( EEG ,'type','modelled','addResiduals',1,'channel',1,cfg{:}); 
        % set(v, 'Position',[   616   707   806   310])
        saveas(gcf,fullfile(dir_pre,unfdir,[pID '_deconv_trialwise.png']) )

        set(0,'CurrentFigure',v); clf
        uf_erpimage( EEG ,'type','residual','channel',1,cfg{:}); 
        saveas(gcf,fullfile(dir_pre,unfdir,[pID '_resid_trialwise.png']) )

        %% plot estimated ERP compare to no-deconvolution approach
        % without deconv
        set(0,'CurrentFigure',v); clf
        g = uf_plotParam(uf_nodc,'channel',1,'deconv',0,'baseline',[-0.1 0],'add_intercept' ,1,'figure',0,'sameyaxis','independent');

        % set nice colors
        for gg = 1:length(g.results.geom_line_handle)
            g.results.geom_line_handle(gg).EdgeColor = [200 0 0]/255;
        end

        % with deconv
        g = uf_plotParam(uf_dc,'channel',1,'deconv',1,'baseline',[-0.1 0],'add_intercept' ,1,'gramm',g,'figure',0,'sameyaxis','independent');
        % set nice colors
        for gg = 1:length(g.results.geom_line_handle)
            g.results.geom_line_handle(gg).EdgeColor = [0 200 0]/255;
        end
        set(v, 'Position',[          31         361        1197         277])

        legend('no deconv','deconv')
        saveas(v,fullfile(dir_pre,unfdir,[pID '_fixR.png']) )
    
        %% save processed EEG
        pop_saveset(EEG, 'filename',pID, 'filepath',fullfile(dir_pre,unfdir), 'savemode','onefile')
        % pop_saveset(dc_FRP, 'filename',[pID '_dc_epochs'], 'filepath',fullfile(dir_pre,unfdir), 'savemode','onefile')

end
