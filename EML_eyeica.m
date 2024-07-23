% Load EEG and align it to eyetracker data
% perform eyetracker-optimised ICA (OPTICAT)
%------------------------------------------

% eyetracker timestamps are milliseconds measured from the time when the tracker software was started
% eeg timestamps are sample number but the sampling rate is 1000Hz here so
% no conversion necessary. **warning** IF reusing this script, you will need to be
% cautious about wehther you are using time in sample or in milliseconds
% and convert where appropriate.

% Load preprocessed data from MNE-python
% - this data has been cleaned using robust detrending, DSS_based line
%   noise removal, and has events labelled from the trial log file.
% - ** It might be worth porting all the MNE-python prepro back to Matlab
%   seeing as it mostly relies on libraries adapted from Matlab in the first
%   place **

clear all; close all
eeglab nogui % sets path defaults
% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
sublist = [67];%TODO random subjects fail e.g. 99
exclude = [20:26 27 31 39 40 78 77 88 138 121 129 160 166]; % Subj to exclude because no eeg or no trigger etc. 54 i fixable - check message file 6. 121,166 not unique?? 
sublist = sublist(~ismember(sublist,exclude) );% & ismember(sublist,find(hasTriggerList.sdcard==1)));

dir_pre = fullfile('..','..','Data','EEG_processed') ;
dir_in = fullfile(dir_pre, 'preprocessed_set');
opticat_dir = 'opticat_cleaned_stricter';
mkdir( dir_pre, opticat_dir)

h_compfig=figure('name','components');
h1=figure('name','Saccade-related potentials before and after OPTICAT-ICA');
h2=figure('name','Fixation-related potentials before and after OPTICAT-ICA');
h3=figure('name','Blink-related potentials before and after OPTICAT-ICA');


for s = 1:length(sublist)
    
    tic;
    pID = ['EML1_',sprintf('%03d',sublist(s))];
    
    %     % use Fieldtrip to load fif data from MNE (TODO can EEGlab read fif
    %     % directly??)
    %     cfg = [];
    %     cfg.dataset = fullfile(dir_fif, [pID '_p.fif']);
    %     EEG_ft = ft_preprocessing(cfg);
    
    EEG = pop_loadset([pID '.set'], dir_in);
    
    % remove fixations and saccades with abnormal params (ask Megan/candace
    % for current cutoffs) In future read from some common file for all of us
    fix_min_dur = 50 /1000*EEG.srate;
    fix_max_dur = 5000/1000*EEG.srate;
    sac_max_dur = 300/1000*EEG.srate;
    fix_keep = strcmp('fixation_either_eye',{EEG.event.type}) & [EEG.event.duration] <= fix_max_dur & [EEG.event.duration] >=fix_min_dur;
    sac_keep = strcmp('saccade_either_eye',{EEG.event.type}) & [EEG.event.duration] <= sac_max_dur ;
    blink_keep = strcmp('blink_either_eye',{EEG.event.type}) ;
    
    event_new = EEG.event(fix_keep | sac_keep | blink_keep);
    EEG.event = event_new;
    
    
    %% OPTICAT ICA
    % This script implements the procedures from:
    % Dimigen, O. (2020). Optimizing the ICA-based removal of ocular artifacts
    % from free viewing EEG. NeuroImage, https://doi.org/10.1016/j.neuroimage.2019.116117
    
    % PARAMETERS:
    opticat_hp = 2; % highpass filter in Hz used for ICA training
    opticat_lp = 100;
    opticat_downsample = 200; % downsample to rate (Hz)
    opticat_epoch_dur = 3; % arbitrary epoch length (s) to cut the data into (because we don't have a straightforward event-related dataset)
    opticat_reject_uV = 300; % reject pseudoepochs from training data with range exceeding this threshold
    opticat_overweighting = 1; % additional spike potentials as proportion of orig data
    opticat_mean_centre = 'true';
    eeg_channels = 1:8; % use all 8 EEG channels - the other 3 are accelerometer
    
    %% Preprocess data
    % filter data 
    EEG = pop_eegfiltnew(EEG,opticat_hp,opticat_lp);
    % resample data 
    EEG = pop_resample(EEG, opticat_downsample);
    % check event structure
    EEG = eeg_checkset(EEG);
    % note: you can apply these just to the training data instead if
    % desired, and apply the ICA filters to the original data, but I think
    % downsampling and filtering the main data here is a good idea
    EEG_copy = EEG; % make a copy of the pre-cleaned data for comparison
    EEG_training = EEG;    
 
    %% Make pseudoepoch event structure
    n_pseudo = floor(EEG.pnts/(opticat_epoch_dur*EEG.srate));
    
    event_pseudo = table();
    event_pseudo.latency =  [0:n_pseudo-1]'*opticat_epoch_dur*EEG.srate+1;
    event_pseudo.type =repmat({'pseudo_epoch'}, height(event_pseudo),1);
    event_pseudo.duration = repmat(opticat_epoch_dur*EEG.srate, height(event_pseudo),1);
    
    event_training = [struct2table(EEG.event); event_pseudo];
    event_training = sortrows(event_training, {'latency'});
    EEG_training.event = table2struct(event_training);
   
    

    
    
    
   
    %% apply pseudoepoching
    EEG_training = pop_epoch(EEG_training, {'pseudo_epoch'}, [0 opticat_epoch_dur]);
    
    %% Reject +- 500 uV pseudoepochs (300 for stricter)
    EEG_training = pop_eegthresh(EEG_training, 1, eeg_channels, -opticat_reject_uV, opticat_reject_uV, 0, opticat_epoch_dur, 1,0);
    % eye events in the event field are already numbered with corresponding
    % pseudoepoch number, which means corresponding eye
    % data will also be removed from training set
    % TODO count trials before and after
    reject_events = ismember([EEG_training.event(:).epoch], find(EEG_training.reject.rejthresh)) ;
    EEG_training = pop_rejepoch(EEG_training,find(EEG_training.reject.rejthresh) ,0);
    
    
    
    %% Overweight spike potentials
    %   Repeatedly append intervals around saccade onsets (-20 to +10 ms)
    %   to training data (-40 to 40 for stricter)
    EEG_training = pop_overweightevents(EEG_training,'saccade_either_eye',[-0.04 0.04],opticat_overweighting,opticat_mean_centre);
    
    %% Run ICA
    fprintf('\nRunning ICA on optimized training data...')
    EEG_training = pop_runica(EEG_training,'chanind',eeg_channels, 'icatype','fastica'); % or use binary ICA for more speed
    
    % Remember ICA weights & sphering matrix
    wts = EEG_training.icaweights;
    sph = EEG_training.icasphere;
    
    %% Transfer unmixing weights
    % Remove any existing ICA solutions from your original dataset
    EEG.icaact      = [];
    EEG.icasphere   = [];
    EEG.icaweights  = [];
    EEG.icachansind = [];
    EEG.icawinv     = [];
    
    fprintf('\nTransfering ICA weights from training data to original data...')
    EEG.icasphere   = sph;
    EEG.icaweights  = wts;
    EEG.icachansind = eeg_channels;
    EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv
    
    %% Eye-tracker-guided selection of ICs
    fprintf('\nIdentifying ocular ICs via saccade/fixation variance-ratio threshold...')
    
    IC_THRESHOLD     = 1.1;   % variance ratio threshold (1.1 was used in in Dimigen, 2020 )
    SACC_WINDOW      = [5 10]; % saccade window (in samples!) to compute variance ratios - before, affter sacc onset
    PLOTFIG          = true;  % plot a figure visualizing influence of threshold setting?
    ICPLOTMODE       = 0;     % plot component topographies (inverse weights)? (2 = only plot "bad" ocular ICs)
    FLAGMODE         = 3;     % overwrite existing rejection flags? (3 = yes)
    MAXCOMP = 3; % how many components to remove MAXIMUM? (3 for stricter)
    
    % Automatically flag ocular ICs (Plöchl et al., 2012)
    [EEG, varratiotable] = pop_eyetrackerica(EEG,'saccade_either_eye','fixation_either_eye',SACC_WINDOW,IC_THRESHOLD,FLAGMODE,PLOTFIG,ICPLOTMODE);
    h=gcf();
    saveas(h, fullfile(dir_pre,opticat_dir, [pID, '_varratios.png']))
    
    % Remove flagged ocular ICs
    badcomps = EEG.reject.gcompreject;
    % choose max 2 with highest variance
    [vr,worst2comps] = maxk(varratiotable(:,3),MAXCOMP);
    badcomps = intersect(find(badcomps), worst2comps);
    
    % plot chosen components & remaining components
    change_current_figure(h_compfig)
    clf
    for ci=1:length(badcomps)
        subplot(1,MAXCOMP,ci)
        topoplot(EEG.icawinv(:,badcomps(ci)),EEG.chanlocs, 'verbose', ...
            'off', 'chaninfo', EEG.chaninfo, 'numcontour', 8);
        title(['Component ' num2str(badcomps(ci)) '| Var ratio = ' num2str(varratiotable(badcomps(ci),3))])
    end
    sgtitle([ pID ' eyeca removed components'],'Interpreter','none')
    saveas(h_compfig, fullfile(dir_pre, opticat_dir, [pID, '_detected_components.png']))
    
    %% remove bad components
    EEG      = pop_subcomp(EEG,badcomps); % remove them
    
    %% view EEG data with ET events
    eegplot(EEG.data,'events',EEG.event, 'data2',EEG_copy.data)
    
    %%  Plot SRP + FRP before and after artefact removal
    % Compute saccade-locked ERPs from clean data, then baseline (ugh mixed units - seconds then milliseconds... and this is SLOW)
    epoch_sac_clean = pop_rmbase(pop_epoch(EEG,{'saccade_either_eye'},[-0.2 0.6]),[-200 0]);
    epoch_sac = pop_rmbase(pop_epoch(EEG_copy,{'saccade_either_eye'},[-0.2 0.6]),[-200 0]);
    epoch_fix_clean = pop_rmbase(pop_epoch(EEG,{'fixation_either_eye'},[-0.2 0.6]),[-200 0]);
    epoch_fix = pop_rmbase(pop_epoch(EEG_copy,{'fixation_either_eye'},[-0.2 0.6]),[-200 0]);
    epoch_blink_clean = pop_rmbase(pop_epoch(EEG,{'blink_either_eye'},[-0.2 0.6]),[-200 0]);
    epoch_blink = pop_rmbase(pop_epoch(EEG_copy,{'blink_either_eye'},[-0.2 0.6]),[-200 0]);
    
    change_current_figure(h1); clf; clear ax
    ax(1)=subplot(1,2,1);
    plot(epoch_sac_clean.times, mean(epoch_sac_clean.data(eeg_channels,:,:),3))
    ylabel('Saccade-related ERP')
    xlabel('Time after saccade [ms]')
    title([pID ': orig data' ], 'Interpreter','none')
    ax(2)=subplot(1,2,2);
    plot(epoch_sac.times, mean(epoch_sac.data(eeg_channels,:,:),3))
    ylabel('Saccade-related ERP')
    xlabel('Time after saccade [ms]')
    title(['clean data - ' num2str(length(badcomps)) ' removed components'])
    legend({EEG.chanlocs.labels}, 'Location','northwest')
    set(h1,'Position',[190 340 1200 400])
    sgtitle('Saccade-related potentials before and after OPTICAT-ICA')
    linkaxes(ax);
    saveas(h1, fullfile(dir_pre, opticat_dir, [pID, '_saccades.png']))
    
    % zoom in on saccade
    set(gca,'XLim',[-20 20]);
    sgtitle('Saccade spike potentials before and after OPTICAT-ICA')
    linkaxes(ax);
    saveas(h1, fullfile(dir_pre, opticat_dir, [pID, '_saccades_SP.png']))
    
    change_current_figure(h2);clf; clear ax
    ax(1)=subplot(1,2,1);
    plot(epoch_fix_clean.times, mean(epoch_fix_clean.data(eeg_channels,:,:),3))
    ylabel('Fixation-related ERP')
    xlabel('Time after fixation [ms]')
    title([pID ': orig data' ], 'Interpreter','none')
    ax(2) =subplot(1,2,2);
    plot(epoch_fix.times, mean(epoch_fix.data(eeg_channels,:,:),3))
    ylabel('Fixation-related ERP')
    xlabel('Time after fixation [ms]')
    title(['clean data - ' num2str(length(badcomps)) ' removed components'])
    legend({EEG.chanlocs.labels}, 'Location','northwest')
    set(h2,'Position',[190 340 1200 400])
    sgtitle('Fixation-related potentials before and after OPTICAT-ICA')
    linkaxes(ax);
    saveas(h2, fullfile(dir_pre, opticat_dir, [pID, '_fixations.png']))
    
    % plot blinks
    change_current_figure(h3); clf;clear ax
    ax(1)=subplot(1,2,1);
    plot(epoch_blink_clean.times, mean(epoch_blink_clean.data(eeg_channels,:,:),3))
    ylabel('Blink ERP')
    xlabel('Time after blink [ms]')
    title([pID ': orig data' ], 'Interpreter','none')
    ax(2)=subplot(1,2,2);
    plot(epoch_fix.times, mean(epoch_blink.data(eeg_channels,:,:),3))
    ylabel('Blink ERP')
    xlabel('Time after blink [ms]')
    title(['clean data - ' num2str(length(badcomps)) ' removed components'])
    legend({EEG.chanlocs.labels}, 'Location','northwest')
    set(h3,'Position',[190 340 1200 400])
    sgtitle('Blink potentials before and after OPTICAT-ICA')
    linkaxes(ax);
    saveas(h3, fullfile(dir_pre, opticat_dir, [pID, '_blinks.png']))
    %% save data and plots
    
    % save clean
    pop_saveset(EEG, 'filename',pID, 'filepath',fullfile(dir_pre, opticat_dir), 'savemode','onefile')
    %     % save raw
    %     pop_saveset(EEG_copy, 'filename',pID, 'filepath',fullfile(dir_pre, 'preprocessed_set'), 'savemode','onefile')
    disp(['Done EyeCA for ' pID])
    
    % Finished ICA for 1 subject!
    load gong
    sound(y,Fs)
    
    toc
end