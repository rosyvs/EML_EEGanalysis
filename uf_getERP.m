function EEG_out = uf_getERP(EEG,varargin)
% This function copies the parts of uf_erpImage that are used to get single
% trial deconvolved estimates but actually returns the EEG structure
% RVS
% DESCRIPTION
%
%Arguments:
% Mandatory
%  cfg.channel (integer):   Which channel(s) ?
%
% Specify to-be-plotted data
%  cfg.method(string): default 'deconv'.
%
%        * 'raw'        : Using raw-data
%
%        * 'modelled'   : Using modelled data (y_hat), usually deconv or
%                         no-deconv version (specified in cfg.datafield)
%                         *Note*: If alignto is specified and no further
%                         remove/keep, then only the modelled data relating
%                         to the alignto-event are plotted
%
%        * 'residual'   : Plot residuals (useful to check for non-modelled
%                         stimulus locked activity)
%
%  cfg.datafield (string):   default 'beta'. Any field of "EEG.unfold",
%                            typically beta or beta_nodc - only in custom
%                            cases something else
%
%  cfg.alignto (String/cell of strings): Which event should mark t=0 in the
%                                        erpimage? (don't confuse with
%                                        sort_alignto). Default is all
%                                        events specified in unfold.X
%  cfg.overlap (boolean): default 0. Re-introduce the overlap? This is
%                         useful to check the modelfit, or useful in cases
%                         where cfg.keep/cfg.remove is used.
%
%  cfg.addResiduals (integer 0-2):  default 0.
%                                   If 1: Adds the overlap-including residuals.
%                                         That is: y_cont - X_dc*beta
%                                   If 2: Adds the corresponding residuals,
%                                         this will be identical to "cfg.method='raw'"
%                                         except when cfg.keep/cfg.remove
%                                         is used
%
%  cfg.winrej (2 integer):   default []. rejection matrix used to fit the
%                           model. necessary if you want to remove the same noisy trials you removed
%                           during modelfitting
% Modify y_hat / modelled data
%
%  cfg.keep (cell):   default {}. Requires a: {}
%                           Possibilty to specify what event-predictor combinations to keep.
%                           In case of "event1:~1+condA, event2:~1+condB"
%                           {'condA','2_(Intercept)'}
%                           would plot single trials without
%                           event1:(Intercept) and event2:condB
%
%  cfg.remove (cell):   default {}. Similar to keep, but specify only the
%                           combinations you want to remove.
%
%  cfg.baseline (2 float)
% Sorting the erpimage
%  cfg.sort_by ():   default ''. XXXmeter is not used.
%
%  cfg.sort_alignto (string/cell of strings):   default 'cfg.alignto'. Event(s) the
%                               ERPimage should be aligned to (specifies
%                               where to look for information for the
%                               black-sorting line). By default alignto
%                               aligns to the erpimage-event (cfg.alignto)
%
%  cfg.sort_by (string):   default 'latency'. Looks in EEG.event.(cfg.by) for
%                           the value to sort by. Has to be an eventname of
%                           EEG.event
%
%  cfg.sort_time (2 integer):   default 'whole epoch'. You can subset where to
%                     look for "cfg.sort_alignto" events. Three typical
%                     use-cases:
%                     1) cfg.time[0 1.5], in case the sort_alignto event
%                        of interest occurs after the cfg.alignto event
%                     2) [-0.3 0] in case the event is before (check out cfg.sort_direction in this case)
%                     3) cfg.time=[0,0], in case the erpimage should be sorted by a field in cfg.alignto
%
%  cfg.sort_direction (string):   default 'forward'. Should the first or
%             last event in cfg.sort_time be selected, in case there are multiple
%             ones? Can be 'forward' or 'backward'. Useful if you are looking for the
%             first cfg.sort_alignto event prior to the cfg.alignto event
%
%

%
%Returns:
%   * (optional) EEG: EEG.data includes The erpimage data, ntimes x ntrials
%   * (optional) sort: the sorting index used to sort the erpimage
%


assert(isfield(EEG,'data'),'uf_erpimage needs the EEG file after uf_glmfit, before uf_condense')
% assert(ismatrix(EEG.data)) % data have to be continuous

if ~isfield(EEG.event,'urevent')
    for e = 1:length(EEG.event)
        EEG.event(e).urevent  = e;
    end
    if isfield(EEG,'urevent') && ~isempty(EEG.urevent)
        error('found EEG.urevent, but not EEG.event.urevent, thus the two things were split up. Fix it by removing EEG.urevent (or setting to []) prior to calling this function')
    end
    EEG.urevent = EEG.event;
end
%% parse & check input
[cfg,EEG_epoch,ufresult] = parse_input(EEG,varargin,nargout);
%% calculate using remove/keep which betas to set to 0
keep = which_parameter_to_keep(ufresult,cfg);
%% calculate data according to remove/keep & the full-beta modelled data
[data,data_yhat] = get_data(ufresult,cfg,keep,EEG_epoch);
%% add all together to a EEG structure, also adds the residuals if necessary
[EEG_out] = get_EEG(ufresult,cfg,data,data_yhat,EEG,EEG_epoch);

% %% calculate the stuff for the erpimage (modified by RVS to skip sortvector
EEG_out = get_sortvector(cfg,EEG_out);
% Subtract the mean value during the baseline window from each trial
if ~isempty(cfg.baseline)
    baseline_idx = (EEG_out.times >= cfg.baseline(1)) & (EEG_out.times < cfg.baseline(2));
    disp(['subtracting mean over ' num2str(sum(baseline_idx)) ' baseline samples' ])
    baseline_mean = mean(EEG_out.data(:, baseline_idx, :), 2);
    EEG_out.data = bsxfun(@minus, EEG_out.data, baseline_mean);
    % assert mean in the baseline is 0 for all atrials
    assert(all(all(abs(mean(EEG_out.data(:, baseline_idx, :), 2)) < 1e-3)), 'Baseline correction failed')
end
end

function [cfg,EEG_epoch,ufresult] = parse_input(EEG,input,n_argout_caller)
cfg = finputcheck(input,...
    {...%#'method','string',{'raw','deconv','nodeconv','residuals','fullmodel','deconvDirect'},'deconv'; %
    'type','string',{'raw','modelled','residual'},'modelled'; % whether to use raw data
    'datafield','string',[],'beta'; % which field to plot, popular: beta & beta_nodc
    'overlap','boolean',[],0;   % re-introduce the overlap
    'addResiduals','integer',[0:2],0; % add residuals?
    ...
    ... %in case of "modelled
    'remove','cell',[],{};     % which event(s) predictor(s) pair(s) to remove {'eventA',{'(Intercept)','stimA'}} (or a cell array of such cell arrays)
    'keep','cell',[],{};       % which event(s) predictor(s) pair(s) to keep   {'eventA',{'stimB'}} (or a cell array of such cell arrays)
    ...
    ... % what to cut & sort to
    'alignto','',[],{};    % what should be the event to align the erpimage?
    'sort_alignto','',[],[];     % default same as alignto (i.e. if you want to search for alignmentevents starting from a different event)
    'sort_by','',[],'latency'; % could be in plotting function? What Eventfield to sort by
    'sort_time','real',[],[];  % when to when to look for the sort_align event
    'sort_direction','string',{'forward','backward'},'forward'; % in case of multiple events, take first or last?
    ...
    ... % Plot Options
    'split_by','',[],[];       % make multiple subplots with separate ERPimages
    'winrej','real',[],[];     % remove parts of the continuous data
    'baseline','real',[],[];
    'plot','boolean',[],n_argout_caller==0;     % if called without requesting output,
    'timelimits','real',[],[]; %from when to when
    'channel','integer',1:size(EEG.data,1),[];
    'figure','boolean',[],1;
    'caxis','real',[],[];
    },'mode','error');


if ischar(cfg); error(cfg);end
assert(isempty(cfg.keep)|isempty(cfg.remove),'either keep or remove, but dont specify both')
assert(~isempty(cfg.channel),'please specify a channel (or set of channels) for the ERPimage')

if ischar(cfg.alignto)
    cfg.alignto = {cfg.alignto};
end
if ~isempty(cfg.sort_alignto) && ischar(cfg.sort_alignto)
    cfg.sort_alignto = {cfg.sort_alignto};
end
if isempty(cfg.sort_alignto)
    cfg.sort_alignto = cfg.alignto;
end

assert(iscell(cfg.alignto))
assert(iscell(cfg.sort_alignto))

%%
if ismatrix(EEG.data)
    % we should epoch the data
    [EEG_epoch] = uf_epoch(EEG,struct('winrej',cfg.winrej,'timelimits',EEG.unfold.times([1,end])+[0 diff(EEG.unfold.times([1:2]))]));
    EEG_epoch = uf_glmfit_nodc(EEG_epoch,'channel',cfg.channel);
    ufresult= uf_condense(EEG_epoch);
else
    EEG_epoch = EEG;
    ufresult= uf_condense(EEG_epoch);

end
assert(isfield(ufresult,cfg.datafield),sprintf('"datafield":%s, not found',cfg.datafield))
end

function [keep] = which_parameter_to_keep(ufresult,cfg)
%% Find which parameters to keep in the model
% default is to keep all
keep = 1:length(ufresult.unfold.colnames);


if ~(isempty(cfg.remove) && isempty(cfg.keep))
    % we fill it only with the keep variables
    list = [];

    if ~isempty(cfg.keep)
        fieldcontent = cfg.keep;
    else
        fieldcontent = cfg.remove;
    end
    for instance = fieldcontent
        ix_variable = find(strcmp(instance{1},ufresult.unfold.variablenames));
        assert(length(ix_variable) == 1,'Could not find variablename:%s, did you misspell something?',instance{1})
        list = [list find(ufresult.unfold.cols2variablenames == ix_variable)];

    end


    if ~isempty(cfg.keep)
        % keep the list only ("define" method)
        disp('Keeping the following predictors:')
        fprintf('%s, ',ufresult.unfold.colnames{list})
        fprintf('\n')
        keep = list;
    else
        % remove the list from all ("subtract" method)
        disp('Removing the following predictors:')
        fprintf('%s, ',ufresult.unfold.colnames{list})
        keep(list) = [];
        fprintf('\n')
    end

end

keep = sort(keep);
end

function [data,data_yhat] = get_data(ufresult,cfg,keep,EEG_epoch)
% Function to calculate yhat=X*b with possibility to remove some of the
% betas (set them to 0)
    function data = calculate_yhat(X,ufresult,cfg,keep)
        beta = ufresult.(cfg.datafield)(:,:,keep);
        if size(ufresult.unfold.Xdc,1) == size(X,1)
            dc=1;
        else
            dc=0;
        end
        assert(~all(isnan(beta(:))),'found only nans, did you specify & calculate the correct channel?')
        % beta = squeeze(mean(beta,1));
        if dc == 1
            data=zeros(length(cfg.channel), size(X,1)); % channel, time
        else
            data=zeros(length(cfg.channel), size(beta, 2),size(X,1));% channel, epochtime, trials
        end
        for c = cfg.channel
            beta_c = beta(c,:,:);
            beta_c = squeeze(mean(beta_c,1));
            if size(beta_c,1) == 1
                beta_c = beta_c';
            end
            if dc == 1
                % timeexpanded
                keepX = ismember(ufresult.unfold.Xdc_terms2cols,keep);
                beta_c = beta_c(:);
            else
                %normal designmat
                keepX = keep;
                beta_c = beta_c'; % because of the "nice" behaviour of matlab to turn every Nx1 to a 1xN we have to undo that
            end
            data_c = X(:,keepX)*beta_c;
            data_c = full(data_c);
            if dc==0
                if size(data_c,1)==1
                    data_c = data_c';
                end
                data(c,:,:) = data_c';
            else
                data(c,:) = data_c';
            end

        end
    end

if strcmp(cfg.type,"raw")
    data = EEG_epoch.data(cfg.channel,:,:);
    data_yhat = nan;
else
    if cfg.overlap
        data      = calculate_yhat(ufresult.unfold.Xdc,ufresult,cfg,keep);
        data_yhat = calculate_yhat(ufresult.unfold.Xdc,ufresult,cfg,1:size(ufresult.unfold.X,2));
    else
        data      = calculate_yhat(ufresult.unfold.X,ufresult,cfg,keep);
        data_yhat = calculate_yhat(ufresult.unfold.X,ufresult,cfg,1:size(ufresult.unfold.X,2));
    end

    if cfg.addResiduals == 1
        % in case of addResiduals == 2 the residuals are already calculated
        data_yhat = calculate_yhat(ufresult.unfold.Xdc,ufresult,cfg,1:size(ufresult.unfold.X,2));
    end
    % % Baseline Correction
    % % Subtract the mean value during the baseline window from each trial
    % if ~isempty(cfg.baseline)
    %     baseline_idx = (ufresult.times >= cfg.baseline(1)) & (ufresult.times < cfg.baseline(2));
    %     baseline_mean = mean(data(:, baseline_idx, :), 2);
    %     data = bsxfun(@minus, data, baseline_mean);
    %     % assert mean in the baseline is 0 for all atrials
    %     assert(all(all(abs(mean(data(:, baseline_idx, :), 2)) < 1e-10)), 'Baseline correction failed')
    % end
    % if ~isempty(cfg.baseline)
    %     data= bsxfun(@minus,data ,mean(data(:,(ufresult.times>=cfg.baseline(1))& (ufresult.times<cfg.baseline(2)),:),2));
    % end
end
end

function [EEG_out] = get_EEG(ufresult,cfg,data,data_yhat,EEG,EEG_epoch)
%% Generate a EEG structure with the new data,
EEG_new = eeg_emptyset;
EEG_new.srate = EEG.srate;
EEG_new.unfold = EEG.unfold;
EEG_new.event = EEG.event;

if ~cfg.overlap || strcmp(cfg.type,"raw")
    % in case of no overlap (or raw), the new data are already epoched
    EEG_out = EEG_epoch;
    EEG_out.data = data;
else
    % first generate continuous data
    EEG_new.data = data;
    % now cut the data once more. This reintroduces overlapping potentials
    % - exactly what we want here
    EEG_out = uf_epoch(EEG_new,'timelimits',ufresult.unfold.times([1,end])+[0 diff(ufresult.times([1:2]))],'winrej',cfg.winrej);

end
%%
% Adding back residuals
if cfg.addResiduals~=0 || strcmp(cfg.type,'residual')
    fprintf('adding residuals\n')
    EEG_residuals = EEG_new;
    EEG_modelled = EEG_new;
    % this is the fully modelled data
    EEG_modelled.data  = data_yhat;
    if ismatrix(data_yhat) % if channel x time => continuous EEG
        %in case of overlap we have continuous residuals and have to cut
        %them once more to epochs
        EEG_modelled= uf_epoch(EEG_modelled,'timelimits',ufresult.unfold.times([1,end])+[0 diff(ufresult.times([1:2]))],'winrej',cfg.winrej);
    end

    % calculate residuals as y-yhat
    EEG_residuals.data = EEG_epoch.data(cfg.channel,:,:) - EEG_modelled.data(cfg.channel,:,:); % this is wrong surely, dims dont add up as yaht here is not epoched !


    if cfg.type == "residual"
        % show only residuals
        EEG_out.data = EEG_residuals.data;
    else
        % add residuals back in
        EEG_out.data = EEG_out.data + EEG_residuals.data;
    end
    clear EEG_residuals % free RAM
end
end

function EEG_out = get_sortvector(cfg,EEG_out)
[keep_epoch]= find(~isnan(eeg_getepochevent(EEG_out,cfg.alignto,[0,0],'type'))); % gets list of epoch indices of event type alignto
assert(~isempty(keep_epoch),'Did not find any events to align to in the epochs')
% in very specific cases, it can happen that we have two different
% eventtypes with exactly the same timing. Because we don't want to double
% epochs, we remove one.
urev = eeg_getepochevent(EEG_out,cfg.alignto,[0,0],'urevent'); % return urevent ix from EEG.urevent table with other events NaN-ed out
urev = urev(keep_epoch); % removes the NaN-ed out values
[~, I] = unique(urev, 'first');
ix_remove = 1:length(urev);
ix_remove(I) = [];
if ~isempty(ix_remove)
    warning('Generated two epochs that occured at the same time - will remove one. This should be happening only rarely.')
end
keep_epoch(ix_remove) = [];
EEG_out.data = EEG_out.data(:,:,keep_epoch);
EEG_out.urevent = EEG_out.urevent(keep_epoch); %added by RVS,
fprintf('Aligning to event %s, %i epochs found\n',strjoin(cfg.alignto,':'),length(keep_epoch)) % RVS caught a bug here.
end

% function [] = plot_erpimage(EEG_out,cfg,sort_vector)
% if isempty(cfg.caxis)
%     cfg.caxis = prctile(EEG_out.data(:),[5 95]);
%     cfg.caxis = [-max(abs(cfg.caxis)) max(abs(cfg.caxis))];
% end
% if isempty(cfg.split_by)
%
%
%     erpimage(EEG_out.data,sort_vector,EEG_out.times,'',10,0,'caxis',cfg.caxis);
% else
% %     [sort_vector,~] = eeg_getepochevent(EEG_out,cfg.alignto,[0,0],cfg.split_by); % output in ms
%     evt_tmp = {EEG_out.urevent.(cfg.split_by)};
%     evt = evt_tmp(cellfun(@(x)~isempty(x),evt_tmp));
%     evt = cellfun(@num2str,evt,'UniformOutput',0); % fixed by Nicolas Langer
%
%     splitlevel = unique(evt);
%     n_splits = length(splitlevel);
%     for n = 1:n_splits
%         subplot(n_splits,1,n)
%         evt_ix = strcmp(evt,splitlevel{n});
%         erpimage(EEG_out.data(:,:,evt_ix),sort_vector(evt_ix),EEG_out.times,'',10,0,'caxis',cfg.caxis);
%
%         title(sprintf('split by:%s, level:%s',cfg.split_by,splitlevel{n}))
%     end
% end
% cmap = cbrewer('div','RdBu',256);
% colormap(cmap(end:-1:1,:));
% end
