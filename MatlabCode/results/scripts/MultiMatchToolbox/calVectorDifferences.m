function M = calVectorDifferences(x1,x2,y1,y2)

% Calculate vector differences (location as well?). Each row of M is the
% euclidean distance between (x1(i), y1(i)) and (x2(i), y2(i)).

%for all the saccades
for k = 1:length(x1)
    
    %calculate Euclidean distance between fixation points
    x_diff = abs(x1(k)*ones(1,length(x2)) - x2);
    y_diff = abs(y1(k)*ones(1,length(y2)) - y2);
    M(k,:) = sqrt(x_diff.^2 + y_diff.^2); 
end



