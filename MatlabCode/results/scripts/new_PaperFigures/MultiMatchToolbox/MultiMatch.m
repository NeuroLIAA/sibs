function MultiMatch(screenSize, varargin)
%
%This script requires the bioinformatics toolbox.
%
%Multimatch assumes you have a 2-level between-subjects manipulation and
%you want to compare the scanpaths of participants in these two conditions
%(e.g.  high and low anxiety people in viewing a set of images). Every item
%will have N number of participants that saw it. Compare all N participants
%with each other on the same item. Then do this for every item. This will
%be outputted in the final data file. If there are 2 images and 3
%participants, p1 and p2 will be compared on image X, p1 vs p3 on imageX,
%p2 vs p3 on imageX, and then the same comparisons will be run for imageY
%(but no comparisons will be made with people viewing imageX versus
%imageY).
%
%The following assumptions about the data have been made:
%
% -all entries must be numeric
% -data must be sorted by participant, then by item,and then by fixation number
% -stimulus IDs are numbers from 1 through to N
%
%
%A matrix of data will be outputted, where each row is a comparison of two
%selected participants on a given trial/item. The columns of the matrix
%will be, in order: (1) stimulus name A, (2) stimulus ID A, (3) participantA 
% number, (4) participantB number, (5) pair type (which conditions 
% participants A and B belong to), (6) pair index (a count of the number 
% of comparisons conducted for each pair type for the current stimulus), 
% (7:11) the outputs of the MultiMatch toolbox, which are:
%
%7) euclidean shape difference 
%8) direction difference
%9) length difference                            type "help doComparison"
%10) position/euclidean distance difference          for more info on these
%11) fixation duration difference
%
% followed by (12) is this participantA-stim pair duplicated? (1=yes). If
% full comparisons are specified (see below), then column (13) will specify
% whether the row for this participantB-stim pair is duplicated, followed
% by (14) stimulus name B and (15) stimulus ID B.
%
% If full comparisons are chosen (see below) then you will be given the
% option to save a summary file. This summary will take the average of all
% the comparisons between a given pair of participants (i.e. between
% participant A and participant B across all stimulus items that the both
% of them saw).
%
% The data will be stored as a .mat file with the name specified by the 
% end user at the end of the script. The variables will also be assigned to
% the MATLAB workspace so that they can be accessed outside of this
% function. A .csv file will also be created of the same name. An error .csv
% file will also be created which reports any warnings that were encountered
% during the comparisons. This should always be checked before analysing
% your data.
%
%
%          --          --          --          --          --
%
%MultiMatch can be called in the following ways:
%
%1) >>MultiMatch(screenSize, pNo, cond, stimID, stimName, fixNo, xPos, ...
%      yPos, fixDur, fullComparisons);
%
%  Each item of a column vector should represent a single fixation to an image,
%  and each column is, in order:
%
%  screenSize: dimensions of screen [x y]
%  pNo: an ID value for the participant fixating the image. [Numeric]
%  *cond: a grouping factor for the participant. E.g. 0 could mean this
%    participant comes from a subgroup of participants who are low in
%    anxiety and 1 could mean they come from a subgroup who are high in
%    anxiety. [Numeric]
%  stimID: a unique value for this stimulus for 1:N stimuli. [Numeric]
%  stimName: a stimulus name that must be numeric (need not be rank order).
%    [Numeric]
%  fixNo: the fixation number for this trial (should start at 1 for each
%  new stimulus presented).
%  xPos: the x position of the current fixation. [Numeric]
%  yPos: the y position of the current fixation. [Numeric]
%  fixDur: the duration of the current fixation. [Numeric]
%  fullComparisons: Do you want to run all possible comparisons or restrict
%    to within-item comparisons? 1 = all, 0 = restricted
%  
%  Note, however, that if you call the function without any inputs you will
%  have to manually specify the inputs later (you will not be able to enter
%  a variable in your MATLAB workspace when asked to enter the (x,y)
%  coordinates, for example), so it is highly advised that you call the
%  function and specify the inputs.
%
%  In addition, the size (in pixels) of the screen on which the data were
%  recorded must be given as an [x y] row vector.
%
% *currently a between-subjects condition must be specified. If you do not
% have a between-subjects condition then simply specify all participants as
% being in the same condition.
%
%
%2) >>MultiMatch(screenSize, dataMatrix[, fullComparisons]);
%
%  The first input is the size of the screen that the data were recorded
%  from specified as an [x y] row vector. The data matrix input is not
%  optional in this version of the call. Each row of dataMatrix should
%  represent a single fixation, and the columns should be, from left to
%  right:
%
%  participant number
%  condition*
%  stimulus ID
%  stimulus name
%  fixation number
%  x position
%  y position
%  fixation duration
%  
%
%  All entries should be numeric.
%
% *currently a between-subjects condition must be specified. If you do not
% have a between-subjects condition then simply specify all participants as
% being in the same condition.
%
%
%3) >>MultiMatch(screenSize, dataFileName[, fullComparisons]);
%
%  Rather than having the data preloaded into MATLAB, you can specify the
%  name of a tab-delimited file (e.g. 'data.dat'). For this to work, it
%  must be in same folder as the MultiMatch toolbox or located in a subfolder
%  within this toolbox. The layout of the data should be as specified in
%  call versions #1 and #2 of this function, above. This is perhaps the
%  easiest way to call the function. Note that this version has only been
%  tested using .dat files. Create a .txt file and either Save As... or
%  just change the .txt suffix to .dat. Make sure you remove any headers
%  from the file before calling it
%
%
%4) >>MultiMatch
%
%  This call will require you to enter all the required details manually.
%  You will be prompted for each. You should enter each vector as a row
%  vector (e.g. [1 2 3]) rather than a column vector (e.g. [1; 2; 3]). The
%  reason is that if you are copying and pasting from another program, it
%  is likely easier to create a row than a columm vector. These are
%  converted to column vectors by this script upon entry.
%
%          --          --          --          --          --
%
% <<KNWON LIMITATIONS: DUPLICATE ENTRIES IN INPUTTED DATA>>
% **Note that any duplicate trials will be treated as separate trials. This
% means that your output file will contain some rows where the participants
% and stimulus ID are the same. It will also run comparisons for that trial
% type with every other possible comparison. So if personX saw stimA twice,
% personX-stimA will be compared to all other comparitors specified (either
% within stimuli or across all stimuli) twice (once for each version of the 
% personX-stimA trial). This isn't a problem because the trials are unique 
% (i.e. they actually did see the trial twice during the experiment).
%
% If they didn't see the trial twice, there are errors in the inputted data. 
% This will lead to a lot of spurious data: not only those where it is
% compared to itself, but duplicate rows where the same stim-participant 
% pairing is compared to all other comparitors. To avoid this issue, ensure 
% there is no spurious data in your inputted data set.
%
% The script will try to identify duplicate entries. If the data are
% organised as described at the top of this help comment (under
% assumptions), it should correctly identify duplicates which are co-located
% in the original data file, but note it will miss any duplicates that are
% separated in the input file. If they are not sorted in the way specified 
% in the help above, it may give spurious results (it will most likely miss
% duplicates). In the datafile a column called 'duplicates' will be '1' if
% it is a duplicate of another row.
%
% Because this script is not designed to deal with duplicate entries that
% should not be in the original data set, it is advised that you run the
% RemoveDuplicateTrials function on your data and then use the resulting
% datafile/MATLAB workspace matrix variable as the input to this function.
% See "help RemoveDuplicateTrials" for more info
%
% <<KNOWN LIMITATIONS: FULL COMPARISONS (no of comparisons)>>
% At the moment if you run full comparisons and there are more than 170
% trials (i.e. more than 2,500 comparisons) the data will be written
% periodically to a .csv file. As a result, the number of comparisons made
% for each condition-pairing (i.e. the 5th column of the outputted csv
% file) will contain zeros. The number of comparisons for each condition
% pairing type can be easily calculated using any program that will filter
% data and allow you to then count the number of unfiltered items. In MS
% Excel, you can enter the formula =COUNTIF([range], "HH") to count the
% number of cells in [range] that contain "HH".
%
% Another small issues with the full comparisons if there are many of them: 
% the error file will contain duplicate entries and on some occasions
% may also say "no reported errors" even though there are reported errors.
% That just means that on a particular run of comparing a single
% participant to all other participants there were no errors.
%
% <<KNOWN LIMITATIONS: FULL COMPARISONS (file name)>>
%
% If there are more than 170 trials, as above, then you will be asked early
% on to specify a file name. If you give the same name as a file that
% already exists, it will append data to that file. Do not attempt to
% overwrite data, simply specify a new name for the file.
%
% <<KNOWN LIMTATIONS: TEXTREAD()>>
% The MATLAB textread() function will be removed from future versions of
% MATLAB. It is currently used here inside twoInputs.m, MultiMatch.m and
% potentially in other scripts too. textscan() will replace textread(), and
% in time this script will require that change to be made. See inside
% the readFromHeadedFile() function for an example of textscan being used.
%
%
%
%   MultiMatch Toolbox created by Marcus Nystrom (Marcus.Nystrom@humlab.lu.se)
%   MultiMatch.m function written by Chris Street
%   (c.street@ucl.ac.uk/cstreet1986@gmail.com) on 20/09/12. 
%
%**************************************************************************

%% obtain data

addpath(genpath(cd)); %add all subfolders to MATLAB's path

%delete old files if they exist
if exist('pMerge.csv', 'file')
    prevState = recycle('off'); % turn recycle off to permanently delete files
    while exist('pMerge.csv', 'file') == 1
        delete('pMerge.csv');
    end
    recycle(prevState); %put recycling back to as it was
end

if nargin == 2 %if two inputs
    
    temp = varargin{1}; %grab input to temporary variable
    data = twoInputs(temp);
    
    %get fullComparison value
    fullComparison = input(...
        'Run all comparisons(1) or restricted to within-stimulus(0)?: ');
    
    %if no input given, run full comparisons by default and inform user
    if length(fullComparison) < 1
        fullComparison = 1;
        fprintf('\nRunning all possible comparisons\n\n');
    end
    
elseif nargin == 3
    
    temp = varargin{1}; %grab input to temporary variable
    data = twoInputs(temp);
    fullComparison = varargin{2}; %third input was fullComparison (0 or 1)

        
    
%10 inputs means each column has been specified separately
elseif nargin == 10
    
    
    data = [varargin{1} varargin{2} varargin{3} varargin{4} varargin{5} ...
        varargin{6} varargin{7}, varargin{8}];
    fullComparison = varargin{9};
    
    
elseif nargin == 0
    
    %put message in command window
    fprintf('\nEnter the following as row vectors...\n\n');
    
    %acquire data
    screenSize = input('size of screen [x y] in pixels: ')';
    pNo = input('participant number vector: ')';
    cond = input('condition vector: ')';
    stimNo = input('stimulus ID vector: ');
    stimName = input('stimulus name vector: ');
    fixNo = input('fixation number vector: ')';
    xPos = input('x position vector: ')';
    yPos = input('y position vector: ')';
    fixDur = input('fixation duration vector: ')';
    fullComparison = input('Run all comparisons(1) or restricted to within-stimulus(0)?: ');
    
    data = [pNo cond stimNo stimName fixNo xPos yPos fixDur];
    
else
    error('either 0, 2, 3 or 10 inputs must be specified. Type "help MultiMatch" for more');
    
end


%fullComparison must be either 0 or 1
if fullComparison ~= 0 && fullComparison ~= 1
    error('Last input to MultiMatch should be 0 (restricted)\nor 1 (full comparisons)');
end

%% Warning to end user

fprintf('This script can take a long time to run. Please be patient\n\n');

%% Isolate items
%Currently all the trials for all participants are stacked on top of each
%other. A single trial here means a single item (e.g. a single visual
%display such as a photo that will generate a scanpath).
%
%To make pairing the items easier, a cell matrix will be created. Each cell
%will contain the data for a single trial/item (pno -> fixdur). Each column
%represents a distinct item, such that each entry under a column is a new
%participant looking at the same stimulus. NOTE: this does not mean that
%the same holds true for the rows. The participant in row 20 of stimulus
%column 1 isn't necessarily the same person as the participant in row 20 of
%stimulus column 2. This is unimportant for our analysis, because the unit
%of analysis is each item, not each participant.

[pairs oneFixation oneFixCounter] = isolateTrials(data, fullComparison);

%% Determine pairings
%Here we have chosen, for every item, to compare every participant that saw
%that item. In practical terms, this means that for every column of the
%"pairs" cell matrix, each row must be compared with every other row.
%However, items from one column must never be compared with items in
%another column (i.e. across stimulus comparisons).
%
% sumdata = [stimName stimID pnoA pnoB pairIndex stimIndex rv(1) rv(2) rv(3) 
%               rv(4) rv(5) duplicate?]

[LL LH HH sumdata sumdataString nStimComparisons oneFixation ...
    nCombo participantMerge periodicWriting saveMatFile] = RunComparisonsEachColumn(...
    pairs, screenSize, fullComparison, oneFixation, oneFixCounter);

%only if we've not been writing all along
if periodicWriting == 0
    %now that the number of comparison types have been counted (LL -> HH), put
    %them into the 5th column of the sumdata variable
    for i = 1:size(sumdata, 1)

        %fill in sumdata(b,5): number of paired combinations
        if strcmp(sumdataString{i}, 'LL') %if condition = "LL"
            sumdata(i,5) = LL; %enter the count of the number of such comparisons
        elseif strcmp(sumdataString{i}, 'LH')
            sumdata(i,5) = LH;
        elseif strcmp(sumdataString{i}, 'HL')
            sumdata(i,5) = LH;
        elseif strcmp(sumdataString{i}, 'HH')
            sumdata(i,5) = HH;
        else
            error('condition needs to be either 00, 01, 10 or 11');
        end

        if fullComparison == 0 %if only comparing within stimuli
            %fill in sumdata(b,6): number of combinations for a given stim
            stimNo = sumdata(i,2);
            sumdata(i,6) = nStimComparisons(stimNo);
        end

    end
end


%% save data

%place results in variables that will be available after the
%function has run.
varData = 'listDataNums';
varDataString = 'listDataText';
varErrors = 'StimParticipantErrs';

%assign variable into the specified name variable
assignin('base', varData, sumdata);
assignin('base', varDataString, sumdataString);
assignin('base', varErrors, oneFixation);

%if data has been saved as we went along
if periodicWriting == 1
    fprintf('\nYour data has been saved as we went along. The name of\nyour datafile is %s\n\n', ...
        saveMatFile);
    
else %if not
    
    %get filename
    saveMatFile = input('choose a name for the full data file\n(optional, default used if left blank): ', 's');

    %if no name is given, make one up using the current clock time and save
    %a .mat file anyway
    if length(saveMatFile)<1
        saveMatFile = sprintf('dataList_%s', datestr(now,'HHMM'));
        fprintf('\n(I am going to save your data anyway... it is called %s)\n', saveMatFile);
    end

    %save to a csv file
    writeResults(sumdata, sumdataString, saveMatFile, fullComparison);

    % write oneFixation errors to a csv file
    writeErrors(oneFixation, saveMatFile);
    
    %save the .mat file 
    save(saveMatFile);
end


%% create summary of the data, averaging by participant
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

if fullComparison == 1 && periodicWriting == 0 %if full comparisons requested but not periodically written
    
    %ask user if they want to create the summary file
    sumPs = input('create summary of participant comparisons? Yes(1) or No(0): ');
    
    if sumPs == 1
        
        [sumPData participantPairs conditionPairs] = createParticipantSummaries(sumdata, ...
               sumdataString, participantMerge); %calculate averages
        
        %write data
        writeSummary(sumPData, participantPairs, conditionPairs, saveMatFile);
        
    end %if sumPs == 1
    
    %save the .mat file now that all the data is present 
    save(saveMatFile);

elseif fullComparison == 1 && periodicWriting == 1 %if periodically written
    
    %ask user if they want to create the summary file
    sumPs = input('create summary of participant comparisons? Yes(1) or No(0): ');
    
    if sumPs == 1
        
        [sumdata sumdataString participantMerge] = readFromHeadedFile(...
            saveMatFile);
        
        
        [sumPData participantPairs conditionPairs] = createParticipantSummaries(sumdata, ...
       sumdataString, participantMerge); %calculate averages
        
        %write data
        writeSummary(sumPData, participantPairs, conditionPairs, saveMatFile);
        
        %save the .mat file now that all the data is present 
        save(saveMatFile);
        
    end %if sumPs == 1
    
end

%% finish


%remove uninformative variables
clear StimParticipantErrs directionSim durationSim lengthSim positionSim ...
    scanpath1 scanpath2 vectorSim;

%N comparisons
if periodicWriting == 1 && sumPs == 0 %no of actual comparisons doesn't work for this scenario
    fprintf('\n\n\nFinished. Potential for %i comparisons, although may not have been possible to conduct them all.\nType "citations" for crediting this work.', sum(nCombo));
else
    fprintf('\n\n\nFinished. Potential for %i comparisons, of which %i were possible to make.\nType "citation" for crediting this work.', sum(nCombo), size(sumdata,1));
end

%where data can be found
fprintf(' data stored in MATLAB workspace\n("listDataText" and "listDataNums") and as .mat and .csv files in this directory\n\n');


%remove pMerge datafile from the disk
if exist('pMerge.csv', 'file')
    delete('pMerge.csv');
end

