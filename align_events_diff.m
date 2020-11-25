function [x_keep, y_keep, offset, lag] = align_events_diff(x,y,tol)
% align two timeseries of irregular, discrete events,
% allowing for missing events in x, y or a combination of both

% intended use-case is aligning EEG triggers to a log file where the
% trigger cable was a bit shite and sometimes shorted and sometimes just
% didn't transmit a pulse at all.
% So there is a monotonic function taking index in x to index in y, but not
% nccessarily linear.

% inputs
% x,y are both timeseries with each element representing sample number of a
% discrete event. length(x) doesn't have to equal length(y).
%
% tol is the maximum allowable timing jitter in samples, which gets used multiple times in this algorithm:
% purposes of estimating offset
% - . if all differences between difference vectors exceeed tol, the corresponding x event is not
% used in the lag estimate,
% - . Following alignment by estimated lag, any remaining events from
% either vector which do not have a corresponding event in the aligned
% sequence are removed

% outputs
% x_keep and y_keep are indices yielding subsets of x and y where a matching was found
% between events following optimisation.
% offset gives the estimated lag (in terms of vector index) of y behind x
% lag gives estimated lag (in terms of timestamps) of y behind x

% ---- Rosy Southwell 2020

if size(x,2)>size(x,1)
    x=x';
end
if size(y,2)>size(y,1)
    y=y';
end

diffx = diff(x);
diffy = diff(y);

% for each element in diffx find its distance to diffy
for i=1:length(diffx)
%             diffmatrix(:,i) = diffx(i)-diffy;
[m,j]=min(abs(diffx(i)-diffy));
    diffmatchings(i,1) = i;
    diffmatchings(i,2) = j;
    resid(i) = diffx(i)-diffy(j); % remaining unnacounted time once diff intervals are matched
end
% get offset (in indices) between diffx and matching in diffy for all i in diffx
all_offsets = diffmatchings(:,2)-diffmatchings(:,1);

% what is modal qualifying offset?
offset=mode(all_offsets(abs(resid)<=tol));

% realign events based on estimated offset
    x_al = x;
    y_al = y;
if offset <0
    y_al=[NaN(-offset,1); y];   
elseif offset>0
    x_al=[NaN(offset,1); x];
end

% NaN pad to get x_al and y_al to same length
if length(x_al)>length(y_al)
    y_al=[y_al; NaN(length(x_al)-length(y_al),1)];
end
if length(y_al)>length(x_al)
    x_al=[x_al; NaN(length(y_al)-length(x_al),1)];
end

% estimate lag based on aligned events
all_lags = y_al-x_al;
lag_candidate=all_lags(find(abs(resid)<=tol)+abs(offset)); % + 1 because diff vector is 1 shorter than orig
lag_candidate_r = tol*round(lag_candidate/tol); % round to nearest tol 
lag_r=mode(lag_candidate_r); % find modal rounded value
% get more precise estimate by averaging over all lags within tol of rounded value
lag=mean(lag_candidate(lag_candidate_r==lag_r));

% align timestamps based on lag
xx = x_al;
yy = y_al-lag;

% find matching events in aligned, offset vectors
xx_keep=[]; yy_keep =[];
xxyydiff = yy-xx;
for k=1:length(xx)
    if any(abs(yy-xx(k))<=tol)
        xx_keep = [xx_keep; k];
        [m, yki] = min(abs(yy-xx(k)));
        yy_keep = [yy_keep; yki];
    end
end

%% deal with duplicates
rem=[];
[~, ind]=unique(xx_keep);
dupx_ix = setdiff(1:length(xx_keep),ind); % indices of duplicated

if ~isempty(dupx_ix) % multiple y values found for given x
    dupx_ev = xx_keep(dupx_ix); % what events are duplicated?
    for ev = unique(dupx_ev)'
   ix=find(xx_keep==ev); % indices in xkeep / ykeep that pertain to this duplicate event in xkeep
        % choose closest or earliest in case of tie
        [m, i] = min( abs( yy(yy_keep(ix)) - xx(xx_keep(ix)) ));
        rem = [rem; ix(setdiff(1:length(ix),i))]; % xx_keep/yy_keep indices to remove. Remove all except the chosen one
    end
end

[~, ind]=unique(yy_keep);
dupy_ix = setdiff(1:length(yy_keep),ind); % indices of duplicated

if ~isempty(dupy_ix) % multiple y values found for given x
    dupy_ev = yy_keep(dupy_ix); % what events are duplicated?
    for ev = unique(dupy_ev)'
        ix=find(yy_keep==ev); % indices in xkeep / ykeep that pertain to this duplicate event in xkeep
        % choose closest or earliest in case of tie
        [m, i] = min( abs( yy(yy_keep(ix)) - xx(xx_keep(ix)) ));
        rem = [rem; ix(setdiff(1:length(ix),i))]; % xx_keep/yy_keep indices to remove. Remove all except the chosen one
    end
end

xx_keep = xx_keep(setdiff(1:length(xx_keep),rem));
yy_keep = yy_keep(setdiff(1:length(yy_keep),rem));

%% get keep indices of orig vector (i.e. account for offset)
x_keep=xx_keep;
y_keep=yy_keep;
if offset <0
    y_keep=yy_keep+offset;
elseif offset>0
    x_keep=xx_keep-offset;
end


% xkeep and ykeep should be the same length, indicating 1-to-1 mapping
assert(length(x_keep)==length(y_keep),'Error: xkeep and ykeep are not the same length!')

% re-estimate lag
lag = mean(y(y_keep)-x(x_keep));

% plots
% take original x,y - not aligned - and plot. 
% Y still needs to be corrected for timestamp lag for visualisatino
% purposes
yyy = y-lag; 
figure(99); clf;
plot(x(x_keep), ones(size(x_keep)),'k.')
hold on
plot(x(setdiff(1:end,x_keep)), ones(length(x)-length(x_keep)),'r.','MarkerSize',18)

plot(yyy(y_keep), zeros(size(y_keep)),'k.')
plot(yyy(setdiff(1:end,y_keep)), zeros(length(yyy)-length(y_keep)),'r.','MarkerSize',18)

line([x(x_keep) yyy(y_keep) ]' ,[ ones(size(x_keep)) zeros(size(y_keep))]')
title('events in x (top) matched to events in y (bottom)')

end
