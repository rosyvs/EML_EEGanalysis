function s = optshift(x,y)

s = 0;
for i=1:length(x)
    for j = 1:length(y)
        s = s - (x(i)-y(j));
    end
end
end % function