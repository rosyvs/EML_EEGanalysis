function [x_keep, y_keep, h] = align_events_time_val(x,y,tol)
% align two timeseries using TIMESTAMPS and event VALUES, both are used to find best match.
% TIMESTAMPS can be approximate match
% VALUES have to be exact match

% allowing for missing events in x, y or a combination of both

% intended use-case is aligning integer trigger values to a log file where
% both have similar timestamps based off the same clock system as is the
% case for fnIRS data - which has regular datetime timestamps

% inputs
% x,y are both timeseries tables with 2 columns:
% TIMESTAMP and VALUE
% timestamp MUST be in the same units for x and y, but is arbitrary - can
% be milliseconds since experiment start, seconds since you woke up, or whatever
%
% length(x) doesn't have to equal length(y).
%
% tol is the maximum allowable timing mismatch in terms of the input timestamp unit
% suggest 1 second for EML

% outputs
% x_keep and y_keep are indices yielding subsets of x and y where a matching was found
% between events following optimisation.
% h is figure handle

% ---- Rosy Southwell 2021



assert(height(x)>2,'x must be a table at least height 2')
assert(height(y)>2,'y must be a table at least height 2')
assert(size(x,2)==2,'x must have 2 columns: TIMESTAMP and VAL')
assert(size(y,2)==2,'y must have 2 columns: TIMESTAMP and VAL')

x.Properties.VariableNames = {'TIMESTAMP','VAL'}; % rename for generalizability
y.Properties.VariableNames = {'TIMESTAMP','VAL'}; % rename for generalizability

% find duration of each event
diffx = [diff(x.TIMESTAMP); 0];
diffy = [diff(y.TIMESTAMP); 0];

% find all possible time lags
for i=1:length(x.TIMESTAMP)
    time_diffmatrix(:,i) = x.TIMESTAMP(i)-y.TIMESTAMP;
    
end

%% find all possible matchings of event value (must be exact)
for i = 1:height(x)
    val_matchmatrix(:,i) =(y.VAL == x.VAL(i));
    
end
val_matchmatrix = double(val_matchmatrix);
val_matchmatrix(val_matchmatrix == 0 )= NaN;
candidate_lags = val_matchmatrix .* time_diffmatrix; % event duration differences for all possible matchings

% round values to nearest tol
candidate_lags_rounded = tol * round(candidate_lags/tol);
% get counts for possible time lags
lag_occurrence = sortrows(tabulate(candidate_lags_rounded(:)),2,'descend');
modal_lag = mean(candidate_lags(abs(candidate_lags-lag_occurrence(1))<tol));
% baseline correct timeseries by modal lag (x-y)
x.TIMESTAMP_ALIGNED = x.TIMESTAMP - modal_lag;


% now find best matches by minimum lag using the VAL based system from
% before
for i = 1:height(x)
    val_matches_in_y = find(y.VAL == x.VAL(i));
    time_matches_in_y = find(abs(milliseconds(y.TIMESTAMP - x.TIMESTAMP_ALIGNED(i)))<tol);
    if ~isempty(val_matches_in_y) && ~isempty(time_matches_in_y)
        candidates = intersect(val_matches_in_y,time_matches_in_y );
        if ~isempty(candidates)
            [mnmm iy]=min(abs(y.TIMESTAMP(candidates) - x.TIMESTAMP_ALIGNED(i)));
            best_match_in_y(i) = candidates(iy);
        end
    else
        best_match_in_y(i) = NaN;
    end
end


best_match_in_y(best_match_in_y==0) = NaN;

x_keep = find(~isnan(best_match_in_y));
y_keep = best_match_in_y(~isnan(best_match_in_y));

% matching must be monotonic!!
while any(diff(y_keep)<0)
    [val,rm] = min(diff(y_keep));
    y_keep=y_keep(setdiff(1:length(y_keep),rm+1));
    x_keep=x_keep(setdiff(1:length(x_keep),rm+1));
    
end

% matching must be one-to-one (no duplicates in y)
if length(unique(y_keep)) < length(y_keep)
    taby = tabulate(y_keep);
    dupey = taby(taby(:,2)>1,1); % dupey is the VALUE that is duplicated - its the index into y
    for d = dupey'
        dupe_loc = find(y_keep==d);% index of the duplicated values in y_keep and x_keep
        dupex = x_keep(y_keep==d); % dupex is the multiple  x indices that have been matched to
        [mnm id]=min(abs(y.TIMESTAMP(d) - x.TIMESTAMP_ALIGNED(dupex)));
        % keep y_keep(dupey(id)) and x_keep(dupex(id)) AKA remove
        y_keep = y_keep(setdiff(1:length(y_keep), setdiff(dupe_loc, dupe_loc(id))));
        x_keep = x_keep(setdiff(1:length(x_keep), setdiff(dupe_loc, dupe_loc(id))));
        
    end
end


if sum(~isnan(y_keep)) >0
    
    % plots
    % take original x,y  and plot.
    
    h=figure(); clf;
    plot(x.TIMESTAMP(x_keep), x.VAL(x_keep),'k.')
    hold on
    plot(x.TIMESTAMP(setdiff(1:end,x_keep)), x.VAL(setdiff(1:end,x_keep)),'r.','MarkerSize',18)
    
    plot(y.TIMESTAMP(y_keep), -y.VAL(y_keep),'b.')
    hold on
    plot(y.TIMESTAMP(setdiff(1:end,y_keep)), -y.VAL(setdiff(1:end,y_keep)),'r.','MarkerSize',18)
    
    line([x.TIMESTAMP(x_keep) y.TIMESTAMP(y_keep) ]' ,[  x.VAL(x_keep)  -y.VAL(y_keep)]')
    title('events in x (top) matched to events in y (bottom)')
else
    h=[]
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Version using datetime %%

%
% function [x_keep, y_keep, h] = align_events_time_val(x,y,tol)
% % align two timeseries using TIMESTAMPS and event VALUES, both are used to find best match.
% % TIMESTAMPS can be approximate match
% % VALUES have to be exact match
%
% % allowing for missing events in x, y or a combination of both
%
% % intended use-case is aligning integer trigger values to a log file where
% % both have similar timestamps based off the same clock system as is the
% % case for fnIRS data - which has regular datetime timestamps
%
% % inputs
% % x,y are both timeseries tables with 2 columns:
% % TIMESTAMP and VALUE
% % timestamp MUST be in the same units for x and y, but is arbitrary - can
% % be milliseconds since experiment start, seconds since you woke up, or whatever
% %
% % length(x) doesn't have to equal length(y).
% %
% % tol is the maximum allowable timing mismatch in terms of the input timestamp unit
% % suggest 1 second for EML
%
% % outputs
% % x_keep and y_keep are indices yielding subsets of x and y where a matching was found
% % between events following optimisation.
% % h is figure handle
%
% % ---- Rosy Southwell 2021
%
%
%
% assert(height(x)>2,'x must be a table at least height 2')
% assert(height(y)>2,'y must be a table at least height 2')
% assert(size(x,2)==2,'x must have 2 columns: TIMESTAMP and VAL')
% assert(size(y,2)==2,'y must have 2 columns: TIMESTAMP and VAL')
%
% x.Properties.VariableNames = {'TIMESTAMP','VAL'}; % rename for generalizability
% y.Properties.VariableNames = {'TIMESTAMP','VAL'}; % rename for generalizability
%
%
% for i = 1:height(x)
%     val_matches_in_y = find(y.VAL == x.VAL(i));
%     time_matches_in_y = find(abs(milliseconds(y.TIMESTAMP - x.TIMESTAMP(i)))<tol);
%     if ~isempty(val_matches_in_y) && ~isempty(time_matches_in_y)
%         candidates = intersect(val_matches_in_y,time_matches_in_y );
%         if ~isempty(candidates)
%             [mnmm iy]=min(abs(y.TIMESTAMP(candidates) - x.TIMESTAMP(i)));
%             best_match_in_y(i) = candidates(iy);
%         end
%     else
%         best_match_in_y(i) = NaN;
%     end
% end
% best_match_in_y(best_match_in_y==0) = NaN;
%
% x_keep = find(~isnan(best_match_in_y));
% y_keep = best_match_in_y(~isnan(best_match_in_y));
%
% if sum(~isnan(y_keep)) >0
%
%     % plots
%     % take original x,y  and plot.
%
%     h=figure(); clf;
%     plot(x.TIMESTAMP(x_keep), x.VAL(x_keep),'k.')
%     hold on
%     plot(x.TIMESTAMP(setdiff(1:end,x_keep)), x.VAL(setdiff(1:end,x_keep)),'r.','MarkerSize',18)
%
%     plot(y.TIMESTAMP(y_keep), -y.VAL(y_keep),'b.')
%     hold on
%     plot(y.TIMESTAMP(setdiff(1:end,y_keep)), -y.VAL(setdiff(1:end,y_keep)),'r.','MarkerSize',18)
%
%     line([x.TIMESTAMP(x_keep) y.TIMESTAMP(y_keep) ]' ,[  x.VAL(x_keep)  -y.VAL(y_keep)]')
%     title('events in x (top) matched to events in y (bottom)')
% else
%     h=[]
% end
% end
