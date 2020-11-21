function [x_aligned, y_aligned, missing_x, missing_y] = align_events_robust(x,y)
% align two timeseries of irregular, discrete events based on inter-event
% timings, allowing for missing events in x, y or a combination of both

% intended use-case is aligning EEG triggers to a log file where the
% trigger cable was a bit shite and sometimes shorted and sometimes just
% didn't transmit a pulse at all.
% So there is a monotonic function taking index in x to index in y, but not
% nccessarily linear.

% inputs
% x,y are both timeseries with each element representing sample number of a
% discrete event. length(x) doesn't have to equal length(y)

% outputs
% x_aligned and y_aligned are subsets of x and y where a matching was found
% between events following optimisation.
% indoces of missing events (wrt
% original timeseries) are returned in missing_x and missing_y

if size(x,2)>size(x,1)
    x=x';
end
if size(y,2)>size(y,1)
    y=y';
end

diffx = diff(x);
diffy = diff(y);
figure(99)
subplot(4,1,1)
plot([x x],[0 1])
xlim([min([x;y]) max([x;y])])
subplot(4,1,2)
plot([y y],[0 1])
xlim([min([x;y]) max([x;y])])

% traverse the matrix
done=0; step=0;i=1; j=1;
while ~done
    
    if i==size(diffx,1) || j==size(diffy,1)
        done=1
        print('finished traverse')
        break
    end
    step=step+1;
    path(step,1)=i;
    path(step,2)=j;
    pathcost(step)=diffx(i)-diffy(j);
    
    % find best next step
    costs(1) = diffx(i+1)+diffx(i)-diffy(j); % skip element in x
    costs(2) = diffx(i)-diffy(j)-diffy(j+1); % skip element in y
    costs(3) = diffx(i+1)-diffy(j+1); % continue diagonally
    [min_cost, direction] = min(abs(costs));
    switch(direction)
        case 1
            i=i+1;
        case 2
            j=j+1;
        case 3
            i=i+1; j=j+1;
    end
    
end

end %fn

% function [x_aligned, y_aligned, missing_x, missing_y] = align_events_robust(x,y)
% % align two timeseries of irregular, discrete events based on inter-event
% % timings, allowing for missing events in x, y or a combination of both
%
% % intended use-case is aligning EEG triggers to a log file where the
% % trigger cable was a bit shite and sometimes shorted and sometimes just
% % didn't transmit a pulse at all.
% % So there is a monotonic function taking index in x to index in y, but not
% % nccessarily linear.
%
% % inputs
% % x,y are both timeseries with each element representing sample number of a
% % discrete event. length(x) doesn't have to equal length(y)
%
% % outputs
% % x_aligned and y_aligned are subsets of x and y where a matching was found
% % between events following optimisation.
% % indoces of missing events (wrt
% % original timeseries) are returned in missing_x and missing_y
%
% if size(x,2)>size(x,1)
%     x=x';
% end
% if size(y,2)>size(y,1)
%     y=y';
% end
%
% xdiffs = diff(x);
% ydiffs = diff(y);
% figure(99)
% subplot(4,1,1)
% plot([x x],[0 1])
% xlim([min([x;y]) max([x;y])])
% subplot(4,1,2)
% plot([y y],[0 1])
% xlim([min([x;y]) max([x;y])])
%
%
%
% for i= 1:length(xdiffs)
%     dists(:,i) = ydiffs - xdiffs(i);
% end
%
% % traverse the matrix
% done=0; step=0;i=1; j=1;
% while ~done
%
%     if i==size(dists,1) || j==size(dists ,2)
%         done=1
%         print('finished traverse')
%         break
%     end
%
%     step=step+1;
%     path(step,1)=i;
%     path(step,2)=j;
%     pathcost(step)=dists(i,j);
%
%
%
%     % find best next step
%     costs(1) = dists(i+1,j); % increment i only
%     costs(2) = dists(i,j+1); % increment j only
%     costs(3) = dists(i+1,j+1); % increment i and j
%     [min_cost, direction] = min(abs(costs));
%     switch(direction)
%         case 1
%             i=i+1;
%         case 2
%             j=j+1;
%         case 3
%             i=i+1; j=j+1;
%     end
%
% end
%
% end %fn



