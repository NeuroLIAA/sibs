function [pairs oneFixation oneFixCounter] = isolateTrials(data, fullComparison)
%
% pairs = isolateTrials(data fullComparison);
%
% This function will take a single matrix of data and separate the
% rows out into separate cells of the "pairs" cell matrix output such that
% each cell contains the fixations for a single item/photo/trial by a single
% participant. Each column will be a specific item, and the rows will be
% participants who saw that item. NOTE that the participant in row 2 col 1
% is not necessarily the same participant in row 2 col 2.
% The input data should have columns that are, from left to right:
%   pNo
%   cond
%   stimName
%   stimID
%   fixNo
%   xPos *
%   yPos * 
%   fixDur *
%
% *necessary for MixMatch comparison (doComparison.m)
%
% ** UPDATE 30/09/12 **
% If a full comparison is requested (=1), then put all the trials into 
% separate cells in a single column, so pairs will be a vector of cells.
% If running a restricted comparison (=0), put into the cell matrix
% described above.
%
% ** UPDATE 08/10/12
% Because full comparisons are now possible, we need to grab the stimulus
% and participant number for both participants involved in the comparison.
% Before this wasn't an issue because the stimulus number was always the
% same (i.e. comparisons were only made within stimulus items, not between
% them).

%predefine pairs variable
pairs = {};
oneFixCounter = 1; %count number of single fixations to item
oneFixation = [0 0 0 0 0]; %predefine variable [stimIDa stimIDb pNoA pNoB error]

z = 1; %start a counter. Allows us to assign to the temporary matrx variable

%take the first row of the matrix into a temporary variable
matrx(z,:) = data(1,:);
z = z+1; %increment counter


for i = 2:size(data,1) %for each row except the first
    
    %if this is the last entry
    if i == size(data,1)
            
        matrx(z,:) = data(i,:); %put last row in current matrix
         
        %if same stimulus, same participant and not first fixation, 'no
        %dupe'
        if data(i,3) == data(i-1,3) && data(i,1) == data(i-1,1) && data(i,5) ~= 1
            matrx(1:size(matrx,1),size(matrx,2)+1) = 0; %set last column to "no duplicate"
        else
            matrx(1:size(matrx,1),size(matrx,2)+1) = 1; %set last column to "no duplicate"
        end
        
        if fullComparison == 0 %within-items, need to determine column number
            y = data(i,3); %get stimID value
        else
            y = 1;
        end
        
        x = findEmpty(y, pairs); %find the next empty row of the cell matrix
        
        %write temporary matrix to "pairs" cell (x, y) of pairs matrix
        pairs{x, y} = matrx;
        clear matrx; %clear the temporary matrix
        
        return; %step out of this function and back to its parent function
        
        
    %if the current stim no (data(i,3)) is the same as the last stim no
    elseif data(i,3) == data(i-1,3)
        
        %check it's the same participant and not the first fixation (could
        %be a duplicate run on the same item)
        if data(i,1) == data(i-1,1) && data(i,5) ~= 1
            
            %take the ith row, all the columns, and place into a temporary
            %matrix
            matrx(z,:) = data(i,:);
            z = z+1; %increment to write to next row of matrx
        
        else %if it's the same stimulus but different participant or the same stim/P but duplicate trial
           
            %if it's a duplicate trial
            if data(i,5) == 1
%                 fprintf('\nDuplicate for stim ID %i, participant %i\n', data(i,3), data(i,1));
                oneFixation(oneFixCounter, 1:5) = [data(i,3) data(i,3) data(i,1)...
                    data(i,1) 4]; %[stim#A stim#B(not applicable) p#A p#B error#4]
                oneFixCounter = oneFixCounter + 1;
                %set all rows and last column of the writing trial to "duplicate"
                matrx(1:size(matrx,1),size(matrx,2)+1) = 1;
            end
            
            if fullComparison == 0 %within-items, need to determine column number
                y = data(i,3); %get stimID value
            else
                y = 1;
            end
            
            x = findEmpty(y, pairs); %find the next empty row of the cell matrix
        
            %write temporary matrix to "pairs" cell (x, y) of pairs matrix
            pairs{x, y} = matrx;
            clear matrx; %clear the temporary matrix
            
            z = 1; %reset temporary matrix counter
            matrx(z,:) = data(i,:); %put new participant into temporary counter
            z = z+1; %increment counter
            
        end
        
        
    elseif data(i,3) ~= data(i-1,3) %if a new stimulus
        
        %no check to see if it's the same participant because that's
        %irrelevent; all we need to know is if it's a new trial,
        %and the fact that the stim number differs means it must be a new
        %trial
         
        if fullComparison == 0 %within-items, need to determine column number
            y = data(i-1,3); %get stimID value of previous trial
        else
            y = 1;
        end
        
        x = findEmpty(y, pairs); %find the next empty row of the cell matrix
        
        matrx(1:size(matrx,1),size(matrx,2)+1) = 0; %set last column to "no duplicate"
        
        %write temporary matrix to "pairs" cell
        pairs{x, y} = matrx;
        
        clear matrx; %clear the temporary matrix
        
        z = 1; %reset temporary matrix counter
        matrx(z,:) = data(i,:); %write this new row to new temp matrix.
        z = z+1;
        
    else
        %something is wrong...
        error('Last stim neither is equal to nor different from current stim... HELP!');
        
    end
    
end

    