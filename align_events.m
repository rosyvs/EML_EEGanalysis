function [sequence, discrepancy] =  align_events(x,y,varargin)
% get the value of the best matching between x and y
% x and y are ordered vectors
% optional argument doplot will plot the input vectors and chosen solution

% Problem outline:
% we can remove elements from x and or y
% a matching comprises subsequences (monotonically increaing index) of x and of y
% of equal length

% we reward keeping elements in x and y by adding a constant V to the value
% for each pair kept
% we decrease the reward for each pairing by the absolute value of the
% difference
if isempty(varargin)
    doplot=false;
else
doplot=varargin{1};
end

if size(x,2)>size(x,1)
    x=x';
end
if size(y,2)>size(y,1)
    y=y';
end

V= 2*max([abs(max(y)-min(x)),  abs(max(x)-min(y)) ] );

m=length(x); n=length(y);

orig_x = x; orig_y=y;
% realign x and y to start at 0;
% x=x-mean(x); y=y-mean(y);


F=zeros(1+m,1+n); % F is the lookup table for our value function
sequences = cell(1+m,1+n);
sequences{1,1}=[];

for i=1:m % loop over x
    for j=1:n % loop over y
        % we are now looking for the sequence representing the best matching up to i and j 
        % consider all previous elements of x
        for k=1:i
           s = F(k,j) + V -abs(x(k)-y(j)); % s is candidate score
           if s>F(i+1,j+1) 
               F(i+1,j+1)=s; 
               sequences{i+1,j+1}=[sequences{k,j}; [k,j]]; %append these indices to prior sequence
           end
        end
        
        for k=1:j-1
           s = F(i,k) + V -abs(x(i)-y(k));
           if s>F(i+1,j+1)
               F(i+1,j+1)=s;
               % store where best value came from
               sequences{i+1,j+1}=[sequences{i,k}; [i,k]];%append these indices to prior sequence

           end
        end
        
    end
end
sequence = sequences{m+1,n+1};
xx=orig_x(sequence(:,1));
yy=orig_y(sequence(:,2));
discrepancy = sum(abs(xx-(yy)));

% diagnistic plot
figure(99); clf
subplot(2,1,1)
plot([orig_x orig_x],[0 1],'r')
hold on
plot([xx xx],[0 1],'k')
xlim([min([orig_x]) max([orig_x])])

subplot(2,1,2)
plot([orig_y orig_y],[0 1],'r')
hold on
plot([yy yy],[0 1],'k')
xlim([min([orig_y]) max([orig_y])])

end % function