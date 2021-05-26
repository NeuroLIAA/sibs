function data = twoInputs(temp)
%determines whether the second argument inputted on the call to
%MultiMatch.m is a file name or a data matrix, and then grabs the
%appropriate data into the 'data' variable.


if isnumeric(temp); %if the second input is numeric, it must be a matrix
        
    data = temp; %take argument into variable "matrx"
        
elseif ischar(temp) %if it's a string
    fileName = temp; %must be the name of the data file

    curdir = cd; %get current directory

    if ismac == 1 %if on a mac
        strcat(curdir, '/', fileName);
    else %otherwise
        strcat(curdir, '\', fileName);
    end

    %read data
    [pNo cond stimNo stimName fixNo xPos yPos fixDur] = ...
        textread(fileName, ' %f %f %f %f %f %f %f %f');

    data = [pNo cond stimNo stimName fixNo xPos yPos fixDur];

else
    error('specify either the variable containing the data matrix or the file name (data.dat)')
end
    