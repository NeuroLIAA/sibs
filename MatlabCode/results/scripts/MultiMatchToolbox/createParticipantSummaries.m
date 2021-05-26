function [sumPData participantPairs conditionPairs] = createParticipantSummaries(sumdata, ...
    sumdataString, participantMerge)
        
%output a summary where every row is a stimulus item, and contains the
%following info. NB It will only write this summary if full comparisons have
%been requested (i.e. across as well as within stimulus items):
%
%   -Participant A
%   -Participant B
%   -Condition pairing
%   -Means of the five outputs from MultiMatch algorithm, averaging across
%   all trials where a comparison between participant A and participant B
%   were made
%   
% sumdata columns:
%
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
%
% It is important that participantMerge is a column vector because
% uniqueInd will be a column vector and the two matrices must be of equal
% dimensions for the averaging to work
 
%% align and get useful data

%make sure participantMerge is a column vector
if size(participantMerge,2) > size(participantMerge,1)
    participantMerge = participantMerge';
end

%determine unique participant pairings
[uniqueID,pMergeIdx,uniqueInd] = unique(participantMerge);

%% check for problems
if size(uniqueInd) ~= size(participantMerge) %there's a problem
    save('lastRun'); %save variables
    error('The number of comparisons is not equal to the size of the dataset. Try deleting the file pMerge.csv and running again. Data saved to lastRun.mat');
end

%if sumdata is larger than the number of comparisons that were run
if size(sumdata, 1) > size(uniqueInd, 1)
    save('lastRun'); %save data
    
    %display error
    s1='The number of data items is larger than the number of comparisons run\n';
    s2='. The error is probably because you typed in a filename that already\n';
    s3=' existed. Choose a new file name and run again. Data from this run is\n';
    s4=' stored in lastRun.mat\n\n';
    err=strcat(s1,s2,s3,s4);
    error(err);
end
    
%% run averages

%average euclidean distance for all unique participant pairings
aveEuc = accumarray(uniqueInd,sumdata(:,7))./accumarray(uniqueInd,1);
%average all other MultiMatch outputs
aveDir = accumarray(uniqueInd,sumdata(:,8))./accumarray(uniqueInd,1);
aveLength = accumarray(uniqueInd,sumdata(:,9))./accumarray(uniqueInd,1);
avePos = accumarray(uniqueInd,sumdata(:,10))./accumarray(uniqueInd,1);
aveDur = accumarray(uniqueInd,sumdata(:,11))./accumarray(uniqueInd,1);


%put into single matrix, with each as a column
sumPData = [aveEuc aveDir aveLength avePos aveDur];

%get the participant pair names into a column vector
participantPairs = uniqueID;

%predefine variable
conditionPairs{length(pMergeIdx)} = '';

%get condition pairing into a column vector
for i = 1:length(pMergeIdx)
    conditionPairs{i} = sumdataString{pMergeIdx(i)};
end
