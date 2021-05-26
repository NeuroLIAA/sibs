function writeTrimmedData(pairs, saveMatFile)

%write the trimmed data from RemoveDuplicateTrials to a .csv file and a
%.dat file. Each cell of pairs contains a matrix. This matrix contains
% the columns:
%
% 1.participant number
% 2.condition (manipulation)
% 3.stimulus name
% 4.stimulus ID (must run from 1:N in increments of 1)
% 5.fixation number of this trial
% 6.x position
% 7.y position
% 8.fixation duration
% 9.less robust attempt to identify duplicates
% 10. Most robust attempt to identify duplicates


%% CSV file (with headers)
%open results file
fileName = strcat(saveMatFile, '.csv');
fid = fopen(fileName, 'w'); %create data file


%create headers
G={'Participant ID', 'Condition', 'Stimulus Name', 'Stimulus ID', ...
    'Fixation Number', 'X position', 'Y position', 'Fixation duration'}; %create headers


%write headers to file
for i = 1:length(G)
    fprintf(fid, '%s,', G{i});
end

fprintf(fid, '\n'); %new line


%get a matrix from the pairs variable
for x = 1:size(pairs,1)
    
    sumdata = pairs{x,1}; %grab the xth matrix
    
    %if 10th column indicates this is not a dupe, write it to the file
    if sumdata(1,10) == 0 
        %for each data entry, write to the file
        for j = 1:size(sumdata, 1) %for no rows

            %write to data file
            fprintf(fid, '%i,', sumdata(j,1:5)); %print relevant integer rows
            fprintf(fid, '%3.2f,', sumdata(j,6:8)); %print relevant float rows
            fprintf(fid, '\n'); %new line

        end %for j =
    end %if
    
end %for x =

fclose(fid); %close the data file

%% .DAT file (without headers, ready for MATLAB import)
%unlike the csv file, separate the columns by tabs rather than commas

%open results file
fileName = strcat(saveMatFile, '.dat');
fid = fopen(fileName, 'w'); %create data file



%get a matrix from the pairs variable
for x = 1:size(pairs,1)
    
    sumdata = pairs{x,1}; %grab the xth matrix
    
    %if 10th column indicates this is not a dupe, write it to the file
    if sumdata(1,10) == 0 
        %for each data entry, write to the file
        for j = 1:size(sumdata, 1) %for no rows

            %write to data file
            fprintf(fid, '%i\t', sumdata(j,1:5)); %print relevant integer rows
            fprintf(fid, '%3.2f\t', sumdata(j,6:8)); %print relevant float rows
            fprintf(fid, '\n'); %new line

        end %for j =
    end %if
    
end %for x =

fclose(fid); %close the data file