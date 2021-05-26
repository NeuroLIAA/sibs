function writeResults(sumdata, sumdataString, saveMatFile, fullComparison, ...
    periodicWriting)
%
%write the sumdata matrix, the number of combinations for each stimulus,
%and the pair index (i.e. number of comparisons for each condition) to a
%csv file.
%
% sumdata columns:
% (1) stimName
% (2) stimID
% (3) participant A
% (4) participant B
% (5) pair index
% (6) stimulus index
% (7) euclidean distance difference
% (8) direction difference
% (9) length difference
% (10) position difference
% (11) duration difference
% (12) duplicate participantA?
%
%
% If full comparisons were made (between as well as within stimuli), we
% need to specify the name and ID of both stimuli in the comparison. In
% this case, sumdata will have 15 columns, with the next two columns being
%
% (13) duplicate participantB?
% (14) stimName B
% (15) stimID B
%
% sumdataString{i} = condition pairing

%if periodic writing isn't defined, it didn't happen
if nargin < 5
    periodicWriting = 0;
end

%open results file
fileName = strcat(saveMatFile, '.csv');

if ~exist(fileName, 'file') %if the file doesn't exist
    
    fid = fopen(fileName, 'w'); %create data file

    if fullComparison == 0
        %create headers
        G={'Stimulus Name', 'Stimulus ID', 'Participant A', 'Participant B', ...
            'Condition Pairing', ...
            'Number of comparisons of type "condition"', ...
            'Number of comparisons for this stimulus', ...
            'Shape difference', 'Direction difference', 'Length difference', ...
            'Position difference', 'Fixation duration difference', 'Duplicate participantA-stim pairing? 1=yes'}; %create headers

    else

        G={'Stimulus Name A', 'Stimulus ID A', 'Participant A', 'Participant B', ...
            'Condition Pairing', ...
            'Number of comparisons of type "condition"', ...
            'Shape difference', 'Direction difference', 'Length direction', ...
            'Position difference', 'Fixation duration difference', ...
            'Duplicate participantA-stim pairing? 1=yes', ...
            'Duplicate participantB-stim pairing? 1=yes', 'Stimulus Name B', ...
            'Stimulus ID B'}; %create headers
    end

    %write headers to file
    for i = 1:length(G)
        fprintf(fid, '%s,', G{i});
    end

    fprintf(fid, '\n'); %new line
    
else %if file already exists
    
    fid = fopen(fileName, 'a+'); %append to data file
    
    %if we shouldn't be appending data due to periodic writing, send a
    %warning
    if periodicWriting == 0
        fprintf('\n\nWARNING: File exists. Appending data to file\n\n');
    end
    
end

    
    
%for each data entry, write to the file
for j = 1:size(sumdata, 1) %for no rows
    
    %write to data file
    fprintf(fid, '%i,', sumdata(j,1:4));
    fprintf(fid, '%s,', sumdataString{j}); %string
    
    if periodicWriting ~= 1

        %if there are full comparisons within and between items, then the 6th
        %column of sumdata will be blank. That's because the number of
        %comparisons conducted for stimulus A was not calculated. This is
        %something that we might want to improve in future versions of MultiMatch
        if fullComparison == 0
            fprintf(fid, '%i,', sumdata(j,5:6)); %integer
        else
            fprintf(fid, '%i,', sumdata(j,5)); %integer
        end
        
    else %periodic writing has been requested
        
        %A downside of running the full comparisons with periodic writing
        %is that it won't tell you how many comparisons of the
        %condition-pairing type have been run. We're just going to fill in
        %the column with 'unknown'
        fprintf(fid, '%s,', 'unknown'); %integer
    end
    
    %write remaining columns (either 7:12 or 7:15)
    fprintf(fid, '%1.3f,', sumdata(j,7:size(sumdata,2))); %3 decimal places
    fprintf(fid, '\n'); %new line
    
end


fclose(fid); %close the data file

