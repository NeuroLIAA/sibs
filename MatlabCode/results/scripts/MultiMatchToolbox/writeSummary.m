function writeSummary(sumPData, participantPairs, conditionPairs, saveMatFile)
% sumPData is a matrix containing the averages for each unique participant
% pairing. The 'averages' are for the 5 output of the MultiMatch script,
% namely (1)shape (2)direction (3)length (4)position and (5)duration.
% Because the participant pair names were turned into a string, they are
% stored separately in participantsPairs{} cell, but the rows will
% correspond to the rows of them matrix
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
% (13) duplicate participantB?
% (14) stimName B
% (15) stimID B
%
% sumdataString{i} = condition pairing

%open results file
fileName = strcat(saveMatFile, '_summary.csv');
fid = fopen(fileName, 'w'); %create data file


G={'[participantA]_[participantB]', 'Condition Pairing', 'Average Shape difference', ...
    'Average Direction difference', 'Average Length difference', 'Average Position difference', ...
    'Average Fixation duration difference'}; %create headers


%write headers to file
for i = 1:length(G)
    fprintf(fid, '%s,', G{i});
end

fprintf(fid, '\n'); %new line


%for each data entry, write to the file
for j = 1:size(sumPData, 1) %for no rows
    
    %write to data file
    fprintf(fid, '%s,', participantPairs{j}, conditionPairs{j}); %string
    fprintf(fid, '%i,', sumPData(j,1:5));
    fprintf(fid, '\n'); %new line
    
end


fclose(fid); %close the data file
