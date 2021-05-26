function x = findEmpty(y, matrx)
%goes through a cell matrix and looks for the next empty row. Unfortunately
%size and length functions don't work on the cell array because even though
%they [i.e. (3,1), (2,2) and (3,2)] are blank, they still exist in the
%following cell matrix:
%
% { [2 x 5]  [2 x 5]  [2 x 5] }
% { [2 x 5]     []    [2 x 5] }
% {   []        []    [2 x 5] }
%
% this function will return, for the above cell matrix:
% x = 3 when y = 1
% x = 2 when y = 2
% x = 4 when y = 3
%
% note that length(matrx(i, :) would always return x = 4
% irrespective of which ith column was selected.

if y > size(matrx,2) %if the column doesn't exist yet
    x = 1; %first row must be blank
    return;
end


for p = 1:size(matrx,1) %for all the rows
    
    temp = matrx{p,y}; %get current row of the specified y column
    
    if length(temp) < 1 %if there's nothing in this cell
        x = p; %p is the xth row that we're looking for
        return; %return to previous function
    end
    
end

%no empty cells, so we need to create a new row. This line won't be reached
%if an empty cell is found in the above loop.
x = size(matrx,1) + 1;
