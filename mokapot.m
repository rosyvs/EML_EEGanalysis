function [sequence, score] =  mokapot(x,y,varargin)
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

m=length(x); n=length(y);
drop_penalty = 500;

E=zeros(1+m,1+n); % lookup table for our 'means' function
V=zeros(1+m,1+n); % lookup table for our 'variance' function
S=zeros(1+m,1+n);
sequences = cell(1+m,1+n);

for i=0:m
    S(1+i,1) = i*drop_penalty;
end
for j=0:n
    S(1,1+j) = j*drop_penalty;
end

for i=1:m % loop over x
    for j=1:n % loop over y
        
        % option 1 is that y(j) is dropped
        newV = V(1+i,j);
        newE = E(1+i,j);
        newS = S(1+i,j) + drop_penalty;
        newsequence = sequences{1+i,j};
           
        % option 2 is that y(j) is matched to x(k) for some k in 1:i
        for k=1:i
%             if k == 1
%                 canV = (x(k)-y(j))^2;
%                 canE = x(k)-y(j);
%                 canS = (i+j-2)*drop_penalty; % s is candidate score
%                 if canS < newS
%                     newV = canV;
%                     newE = canE;
%                     newS = canS;
%                     sequences{i,j}=[1,j]; % append these indices to prior sequence
%                 end
%             else
                canV = V(k,j) + (x(k)-y(j))^2;
                canE = E(k,j) + x(k)-y(j);
                canS = S(k,j) - 2*(x(k)-y(j))*E(k,j) + (i-k)*drop_penalty; % s is candidate score
                if canS < newS
                    newV = canV;
                    newE = canE;
                    newS = canS;
                    newsequence=[sequences{k,j}; [k,j]]; % append these indices to prior sequence
                end
%             end
        end
        
        E(1+i,1+j) = newE;
        V(1+i,1+j) = newV;
        S(1+i,1+j) = newS;
        sequences{1+i,1+j}= newsequence;
    end
end
sequence = sequences{end};
score = S(end);

end % function