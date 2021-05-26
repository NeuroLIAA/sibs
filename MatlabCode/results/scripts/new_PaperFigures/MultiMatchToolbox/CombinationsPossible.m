function [sumdata sumdataString LL LH HH stimCompCounter oneFixation, ...
    oneFixCounter b participantMerge periodicWriting saveMatFile] = ...
    CombinationsPossible(b, oneFixCounter, oneFixation, LL, LH, ...
    HH, p, pairs, N, stimCompCounter, tempx, screenSize, sumdata, ...
    sumdataString, fullComparison, participantMerge, periodicWriting)
%
%This function gets called when there are enough rows in a given column of
%'pairs' to start making comparisons between them. Grab two cells from the
%pairs variable. This function will terminate when all comparisons from a
%single column have been performed
%
%if periodicWriting == 1 then write to a .csv file after comparator A has
%been compared to all other participants.
        
%% predefine variables
breakVal = 0; 


%specify file name if we're going to periodically write to it
if periodicWriting == 1
    %get filename
    saveMatFile = input('choose a name for the full data file\n(optional, default used if left blank): ', 's');

    %after every comparison with Comparator A has been completed, a bit of
    %info will be displayed if set to 1
    verbse = 1;
    
    %if no name is given, make one up using the current clock time and save
    %a .mat file anyway
    if length(saveMatFile)<1
        saveMatFile = sprintf('dataList_%s', datestr(now,'HHMM'));
        fprintf('\n(I am going to save your data anyway... it is called %s)\n', saveMatFile);
    end
    
else
    saveMatFile = '';
    %after every comparison with Comparator A has been completed, a bit of
    %info will be displayed if set to 1
    verbse = 0; 
end

%check the file doesn't already exist. If it does, give a warning
if exist(strcat(saveMatFile, '.csv'), 'file') == 1
    fprintf('\n\nWARNING: File already exists. Appending data to that file. Ctrl + C to stop.\n\n');
end

%% run comparisons

for j = 1:N-1 %for all potential participant xs (comparator As)

%     %DEBUGGING - tells you how many comparisons have been made so far
%     disp(b);
    
    %update x to next value
    x = tempx + 1; %start x from new position
    tempx = x; %reset starting position of x

    for q = x+1:N % all potential participant x's (comparator A) to be paired with all other 
                   %participants (subtract any Ps with only 1 fixation)

        %rather than using q to go from 2 up to N, convert it so we can go
        %from N down to 2. More generally, convert so we can go from N down
        %to x+1. This line looks a bit crazy, but it works...
        Q = (N+2)-q + (j-1);

        if Q > x %only run comparisons if Q (comparator B) is greater than x (comparator A)
            
            %*******************
            %   COMPARATOR A
            %*******************
            data1 = pairs{x, p}; %column p, participant x

            %until data1 has more than 1 fixation
            if size(data1,1) == 1

                %make record [stimID stimIDb pnoA pnoB error1];
                oneFixation(oneFixCounter, 1:5) = [data1(1,3) 0 data1(1,1) 0 1];
                oneFixCounter = oneFixCounter+1;

                %this variable will skip the rest of this loop and restart the
                %'for' loop at "for q = x+1:N"
                breakVal = 1;

            end


            %*******************
            %   COMPARATOR B
            %*******************
            data2 = pairs{Q, p}; %column p, participant Q

            %check that data2 has more than one fixation point (make
            %warning separate from while loop so that the warning is
            %only shown once
            if size(data2,1) == 1
  
                %make record [stimID stimIDb pnoA pnoB error1];
                oneFixation(oneFixCounter, 1:5) = [data2(1,3) 0 data2(1,1) 0 1];
                oneFixCounter = oneFixCounter+1;
            end

            
            %keep updating the data2 variable until it contains more than 1 row. If it
            %has only one row, it means that this particular trial only has one
            %fixation point

            while size(data2,1) == 1
                Q = Q-1; %get next participant value
                if Q <= x+1 %there aren't any more pairs for participant x
                    %make data2 large enough to break the "while" loop
                    data2 = [0 0; 0 0];
                end
                breakVal = 1;
            end



            if breakVal ~= 1 %if not set to 1, then continue

                fixationData1 = data1(:,6:8); %(x, y, dur)
                fixationData2 = data2(:,6:8); %(x, y, dur)


                %********************************
                %         MULTIMATCH
                %********************************

%                %for debug, let's see how many comparisons have been
%                %made so far and between which pairs{stim,participant}
%                fprintf('comparison %i, pairs{%i,%i} vs pairs{%i,%i}\n\n', ...
%                   b,x,p,Q,p);


                suppressOutput = 1; %don't show the results of every call to doComparison

                %run the MultiMatch toolbox
                rv = doComparison(fixationData1, fixationData2, screenSize, ...
                    suppressOutput); 
                %passing the save name creates data-saving
                %redundancy. Can add it as last input to the above
                %function if you wish, but it will delay run time
                %significantly!


                %if only one value in rv, it must be NaN resulting from
                %error inside doComparison. If there are 5 elements to
                %rv, then the doComparison score worked
                if length(rv) == 5 

                    %*******************
                    %    COMPILE DATA
                    %*******************

                    %columns, in order, for final output:
                    sumdata(b,1) = data1(1,4); %stimName
                    sumdata(b,2) = data1(1,3); %stimID
                    sumdata(b,3) = data1(1,1); %participantA
                    sumdata(b,4) = data2(1,1); %participantB

                    %convert the condition value (0 or 1) into strings and then
                    %concatenate


                    if data1(1,2) == 0
                        condA = 'L';
                    elseif data1(1,2) == 1
                        condA = 'H';
                    end

                    if data2(1,2) == 0
                        condB = 'L';
                    elseif data2(1,2) == 1
                        condB = 'H';
                    end

                    condition = strcat(condA, condB); %should be either 00, 01, 10 or 11

                    sumdataString{b} = condition; %write to sumdataString

                    %sumdata(b,5) is defined later because it can only be written when
                    %all the values have been written to the data file
                    %(participantindex)
                    %sumdata(b,6) is defined later for similar reason
                    %(stimindex)
                    sumdata(b,7) = rv(1); %Euclidean distance difference
                    sumdata(b,8) = rv(2); %direction difference
                    sumdata(b,9) = rv(3); %length difference
                    sumdata(b,10) = rv(4); %position difference
                    sumdata(b,11) = rv(5); %duration difference
                    sumdata(b,12) = data1(1,9); %duplicate participant-stim A?
                    
                    if fullComparison == 1 %If full comparisons made, need to grab stim name of second stimulus
                        sumdata(b,13) = data2(1,9); %duplicate participant-stim B?
                        sumdata(b,14) = data2(1,4); %stimName
                        sumdata(b,15) = data2(1,3); %stimID

                        %merge participant numbers, useful if running full
                        %comparisons and a summary file is wanted
                        pA = num2str(data1(1,1));
                        pB = num2str(data2(1,1));
                        participantMerge{b} = strcat(pA,'_',pB);
                        
                    end

                    

                    b = b+1; %increment counter, write the next comparisons to a new line of the output

                    
                    if periodicWriting == 0 %don't bother testing if we're periodically writing
                        if strcmp(condition, 'LL')
                            LL = LL + 1;
                        elseif strcmp(condition, 'LH')
                            LH = LH + 1;
                        elseif strcmp(condition, 'HL')
                            LH = LH + 1;
                        elseif strcmp(condition, 'HH')
                            HH = HH + 1;
                        else
                            error('condition needs to be either 00, 01, 10 or 11');
                        end
                    end
                    
                    %number of stimulus comparisons made
                    stimCompCounter = stimCompCounter + 1;
                        

                else %if length(rv) == NaN

                    %make record [stimID stimIDb pnoA pnoB error2];
                    oneFixation(oneFixCounter, 1:5) = [data1(1,3) data2(1,3)...
                        data1(1,1) data2(1,1) 2];
                    oneFixCounter = oneFixCounter+1;


                end %if length(rv) == 5


            end %if breakVal ~= 1

            breakVal = 0; %reset breakVal

        end %if Q > x 

    end %for q = x+1:N (for all comparisons that can be made with comparator A)
    
    %now that all comparisons for comparator A have been run, if we're
    %writing periodically this is the place to do it. THe data writing
    %technique is to take all the sumdata so far and write it to the file,
    %but if the file already exists then append the data. So for each run,
    %we need to clear sumdata so that we don't append all the old data to
    %the file as duplicates.
    if periodicWriting == 1 && size(sumdata,1) > 1 %if any comparisons were run
        
        
        %save to a csv file
        writeResults(sumdata, sumdataString, saveMatFile, fullComparison, ...
            periodicWriting);
        
        %keep track of a merge of the first participant ID and the second
        %participant ID for each comparison
        saveParticipantMerge(participantMerge);

        % write oneFixation errors to a csv file
        writeErrors(oneFixation, saveMatFile);
        
        %reset variables
        clear sumdata sumdataString oneFixCounter oneFixation b participantMerge;
        sumdata(1,1:12) = 0;
        sumdataString = '';
        oneFixCounter = 1; %count number of errors
        oneFixation = [0 0 0 0 0]; %predefine variable [stimIDa stimIDb pNoA pNoB error]
        b = 1; %start writing into sumdata at b, i.e. at first row
        participantMerge = '';
        
        if verbse == 1
            fprintf('\nNth stimulus-participant pair compared. N = %i', j);
        end
        
    end
    


end %for j = 1:N-1