function pairs = dupeCellDelete(pairs)

%Makes a comparison of each cell with every other cell. If it finds that
%the participant number and stimulus ID in one cell is identical to that of
%another cell, it adds a 1 to the 9th column of the matrix inside the cell.
%Otherwise, it adds a 0. Later this can be used to filter data.
%
%because we're piggy-backing off the isolateTrials function, there will be
%an additional 9th row that attempted to identify duplicates. However, it
%is untested and certainly not as robust as this attempt to remove
%duplicates. So we'll create a 10th row that again attempts to identify
%duplicates in a far more reliable fashion.

x = 0; %start a counter for the second person to pair with participant 
    %q in the loop below
tempx = x; %keep track of the starting position of x
N = size(pairs,1); %number of rows/trials

for j = 1:N-1 %for all potential participant x's

    %update x to next value
    x = tempx + 1; %start x from new position
    tempx = x; %reset starting position of x

    for q = x+1:N % all potential participant x's to be paired with all other participants (subtract any Ps with only 1 fixation)

        %rather than using q to go from 2 up to N, convert it so we can go
        %from N down to 2. More generally, convert so we can go from N down
        %to x+1. This line looks a bit crazy, but it works...
        Q = (N+2)-q + (j-1);

        if Q > x %only run comparisons if Q is greater than x
            
            
            %% get two trials to compare
            %*******************
            %   COMPARATOR A
            %*******************
            data1 = pairs{x, 1}; %column 1, participant x


            %*******************
            %   COMPARATOR B
            %*******************
            data2 = pairs{Q, 1}; %column 1, participant Q

            
            pnoA = data1(1,1);
            pnoB = data2(1,1);
            stimA = data1(1,3);
            stimB = data2(1,3);
            
                
            %% identify duplicates
            
            %duplicate
            if pnoA == pnoB && stimA == stimB %if pnos and stimIDs match
                
                %COMPARATOR A
                if size(data1,2) == 9 %if it hasn't been given a value yet
                    data1(:,10) = 0; %not a duplicate
                end
                
                %COMPARATOR B
                if size(data2,2) == 9 %no value assigned yet
                    data2(:,10) = 1;
                %if it has been identified previously as a non-duplicate
                %but here has been found to be a duplicate, update the
                %value to 1
                elseif data2(1,10) == 0 
                    data2(:,10) = 1; %is a duplicate
                end
                    
                
            else %no duplicate
                
                
                %COMPARATOR A
                if size(data1,2) == 9 %if it hasn't been given a value yet
                    data1(:,10) = 0; %not a duplicate
                end
                
                %COMPARATOR B
                if size(data2,2) == 9 %if it hasn't been given a value yet
                    data2(:,10) = 0; %not a duplicate
                end
                
            end
    
            
            %% write data back to cells
            
            pairs{x, 1} = data1;
            pairs{Q, 1} = data2;


        end %if Q > x 

    end %for q = x+1:N

end %for j = 1:N-1