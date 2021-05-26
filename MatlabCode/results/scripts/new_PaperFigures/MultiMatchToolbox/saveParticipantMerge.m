function saveParticipantMerge(participantMerge)

%% open results file
fileName = 'pMerge.csv';

if ~exist(fileName, 'file') %if the file doesn't exist
    
    fid = fopen(fileName, 'w'); %create data file
    
else %if file already exists
    
    fid = fopen(fileName, 'a+'); %append to data file
    
end
    
    
%% for each data entry, write to the file
for j = 1:size(participantMerge, 2) %for no entries
    
    %write to data file
    fprintf(fid, '%s', participantMerge{j});
    fprintf(fid, '\n'); %new line
    
end

%% close file

fclose(fid); %close the data file