function [x_keep, y_keep, h] = align_events_val(x,y,tol)
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
% length(x) doesn't have to equal length(y).
%
% tol is the maximum allowable timing mismatch in milliseconds,

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


for i = 1:height(x)
    val_matches_in_y = find(y.VAL == x.VAL(i));
    time_matches_in_y = find(abs(milliseconds(y.TIMESTAMP - x.TIMESTAMP(i)))<tol);
    if ~isempty(val_matches_in_y) && ~isempty(time_matches_in_y)
        candidates = intersect(val_matches_in_y,time_matches_in_y );
        if ~isempty(candidates)
            [mnmm iy]=min(abs(y.TIMESTAMP(candidates) - x.TIMESTAMP(i)));
            best_match_in_y(i) = candidates(iy);
        end
    else
        best_match_in_y(i) = NaN;
    end
end
best_match_in_y(best_match_in_y==0) = NaN;

x_keep = find(~isnan(best_match_in_y));
y_keep = best_match_in_y(~isnan(best_match_in_y));

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
