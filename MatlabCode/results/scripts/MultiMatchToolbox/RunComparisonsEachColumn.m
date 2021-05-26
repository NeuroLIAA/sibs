function [LL LH HH sumdata sumdataString nStimComparisons oneFixation ...
    nCombo participantMerge periodicWriting saveMatFile] = RunComparisonsEachColumn(...
    pairs, screenSize, fullComparison, oneFixation, oneFixCounter)
%runs through each column of pairs and makes the multimatch
%comparisons. If full comparisons have been requested, there may be
%hundreds of thousands of even millions of comparisons being requested. If
%there are over 2,500 comparisons (i.e. if MATLAB is unable to calculate
%the factorial values to determine the permutations) then allow for a
%different way to run the script. There are two people in a comparison:
%after comparison the first person to all other participants, write the
%data. Then do the same for the second, and the third, and so on.
%
%MATLAB might run out of memory space or the user may fear a runaway loop,
%hence this periodic writing strategy.

%% predefine variables

b = 1; %initiate counter. Number of rows in final data output.
breakVal = 0; %is set to 1 if a participant has only one fixation point in their scanpath

%count the number of 
LL = 0;
LH = 0;
HH = 0;

sumdata = 0;
sumdataString = '';
participantMerge = '';
periodicWriting = 0;

nCombo = zeros(1, size(pairs,2));
nStimComparisons = zeros(1, size(pairs,2));

%% loop each column
for p = 1:size(pairs,2) %for all the columns of the "pairs" matrix (i.e. for each stimulus)
        
    
    %determine the actual length of the column (i.e. ignore any blank
    %cells) plus 1 (below function actually finds the index for the next
    %blank cell, so we need to subtract 1 from that to know the number of
    %participants). 
    empt = findEmpty(p, pairs);
    N = empt - 1;
    
    if N == 2
        nCombo(p) = 1; %only one viable combination, but factorials will return infinite
    elseif N <= 1
        nCombo(p) = 0;
        oneFixation(oneFixCounter, 1:5) = [p 0 0 0 3]; %[stim# stim#B noParticipant noParticipantB error#3]
        oneFixCounter = oneFixCounter + 1;
    elseif N > 170 %unable to determine the number of combinations; the factorial is just too big
        nCombo(p) = Inf; % 4,988 trials = 12.5 million combinations
        disp('over 2,500 comparisons are being run, please be patient');
        
        %if running full comparisons, allow for period writing of data
        if fullComparison == 1
            periodicWriting = 1;
            
            %delete any pre-existing merge file
            if exist('pMerge.csv', 'file')
                delete('pMerge.csv');
            end
        end
        
    else
        %determine the number of possible combinations of 2 rows (w/o replacement)
        nCombo(p) = factorial(N)/(factorial(2)*factorial(N-2));
                
    end
    
    stimCompCounter = 0; %start counter for number of comparisons run within this stimulus item
    x = 0; %start a counter for the second person to pair with participant 
    %q in the loop below
    tempx = x; %keep track of the starting position of x
 
    
    if nCombo(p) > 0%only allow if there are combinations to run
        
        [sumdata sumdataString LL LH HH stimCompCounter oneFixation, ...
            oneFixCounter b participantMerge periodicWriting saveMatFile...
            ] = CombinationsPossible(b, oneFixCounter, oneFixation, LL, LH, ...
            HH, p, pairs, N, stimCompCounter, tempx, screenSize, sumdata, ...
            sumdataString, fullComparison, participantMerge, periodicWriting);
        

    end %if nCombo > 1
    
    nStimComparisons(p) = stimCompCounter; %number of stimulus comparisons for this item
    
end %for p = 1:size(pairs,2).  i.e. for all stimuli
    
clear empt; %clear temporary variable

