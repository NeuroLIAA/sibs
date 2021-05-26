function writeErrors(errorvals, saveFile)
%
%write the oneFixation matrix (containing errors) to a
%csv file.
%
% oneFixation columns:
%(1) stimulus ID (not name) of comparator A
%(2) stimulus ID (not name) of comparator B
%(3) participant number A
%(4) participant number
%(5) error description (see lines 34-40 below)
%

%% prepare file

if nargin < 2
    saveFile = sprintf('errorList_%s', datestr(now,'HHMM'));
end

%open results file
fileName = strcat(saveFile, '_errors.csv');

%% if file doesnt exist
if ~exist(fileName, 'file') %if the file doesn't exist
    fid = fopen(fileName, 'w'); %create data file

    %headers
    %create headers
    G={'Stimulus ID A', 'Stimulus ID B', 'Participant number A', ...
        'Participant number A', 'Error description'}; %create headers

    %write headers to file
    for i = 1:length(G)
        fprintf(fid, '%s,', G{i});
    end

    
else
    %% if file exists
    fid = fopen(fileName, 'a+'); %append to file
end

fprintf(fid, '\n'); %new line

%% get unique row values

errorValues = unique(errorvals,'rows');


%% write new rows to data file

%for each data entry, write to the file
for j = 1:size(errorValues, 1) %for no rows
    
    
    %write to data file
    fprintf(fid, '%i,', errorValues(j,1:4));
    
    if errorValues(j,5) == 1
        errorDescription = 'only one fixation for this trial';
    elseif errorValues(j,5) == 2
        errorDescription = 'all fixations were out of screen boundaries';
    elseif errorValues(j,5) == 3
        errorDescription = 'not enough participants to make comparisons for this stimulus';
    elseif errorValues(j,5) == 4
        errorDescription = 'Duplicate trial for this participant and stimulus. This participant/trial will be compared to its duplicate in your data file';
    
        %if only one row, which is blank, then there were no errors
    elseif size(errorValues,1) == 1 && errorValues(1,1) == 0 && errorValues(1,2) == 0 ...
            && errorValues(1,3) == 0 && errorValues(1,4) == 0
        errorDescription = 'no errors reported';
    else
        errorDescription = sprintf('unknown, error number %i', errorValues(j,5));
    end
    
    fprintf(fid, '%s,', errorDescription); %string
    fprintf(fid, '\n'); %new line
    
end


fclose(fid); %close the data file

