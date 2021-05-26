function RemoveDuplicateTrials(varargin)
%
%RemoveDuplicateTrials(varargin);
%
% Will take a matrix of data (or a .dat filename that contains a matrix of
% info, tab or comma separated) where each row represents a fixation and
% each column represents, from left to right:
%
% 1.participant number
% 2.condition (manipulation)
% 3.stimulus name
% 4.stimulus ID (must run from 1:N in increments of 1)
% 5.fixation number of this trial
% 6.x position
% 7.y position
% 8.fixation duration
%
% Any single trial is thus represented as a set of rows in this matrix,
% from fixation number (column 5) 1 to N. Duplicate trials will be removed
% from the data file and saved back to a .csv file with headers, a .dat
% file without headers (i.e. data in a format that MATLAB can use in a call
% to MultiMatch) and a .mat file. It is the first entry of the duplicates
% that is kept and other entries that are removed
%
% the following assumptions about the data have been made:
%
% -all entries must be numeric
% -data must be sorted by participant, then by item, and then by fixation number
% -stimulus IDs are numbers from 1 through to N
%
%
%RemoveDuplicateTrials can be called in the following ways:
%
%1) >>RemoveDuplicateTrials(pNo, cond, stimID, stimName, fixNo, xPos, ...
%      yPos, fixDur);
%  Each item of a column vector should represent a single fixation to an image,
%  and each column is, in order:
%
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
%
%
% *currently a between- or within-subjects condition must be specified. If
% you do not have a manipulation then simply specify all participants as
% being in the same condition.
%
%
%2) >>RemoveDuplicateTrials(dataMatrix);
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
%3) >>RemoveDuplicateTrials(dataFileName);
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
%4) >>RemoveDuplicateTrials;
%
%  This call will require you to enter all the required details manually.
%  You will be prompted for each. You should enter each vector as a row
%  vector (e.g. [1 2 3]) rather than a column vector (e.g. [1; 2; 3]). The
%  reason is that if you are copying and pasting from another program, it
%  is likely easier to create a row than a columm vector. These are
%  converted to column vectors by this script upon entry.
%
%  Note, however, that if you call the function without any inputs you will
%  have to manually specify the inputs later (you will not be able to enter
%  a variable in your MATLAB workspace when asked to enter the (x,y)
%  coordinates, for example), so it is highly advised that you call the
%  function and specify the inputs.
%
%   RemoveDuplicateTrials and all subfunctions written by Chris Street
%   (c.street@ucl.ac.uk/cstreet1986@gmail.com) on 30/09/12. 
%
%**************************************************************************


%% prepare to run
addpath(genpath(cd)); %add all subfolders to MATLAB's path

fullComparison = 1; %make 'pairs' a vector, not a matrix. Set to 0 to make it a matrix

disp('This script can take a long time to run with large datasets. Please be patient');

%% obtain data

if nargin == 1 %if two inputs
    
    temp = varargin{1}; %grab input to temporary variable
    data = twoInputs(temp);
        
    
%8 inputs means each column has been specified separately
elseif nargin == 8
    
    
    data = [varargin{1} varargin{2} varargin{3} varargin{4} varargin{5} ...
        varargin{6} varargin{7}, varargin{8}];
    
    
elseif nargin == 0
    
    %put message in command window
    fprintf('\nEnter the following as row vectors...\n\n');
    
    %acquire data
    pNo = input('participant number vector: ')';
    cond = input('condition vector: ')';
    stimNo = input('stimulus ID vector: ');
    stimName = input('stimulus name vector: ');
    fixNo = input('fixation number vector: ')';
    xPos = input('x position vector: ')';
    yPos = input('y position vector: ')';
    fixDur = input('fixation duration vector: ')';
    
    data = [pNo cond stimNo stimName fixNo xPos yPos fixDur];
    
else
    error('either 0, 1, or 8 inputs must be specified. Type "help RemoveDuplicateTrials" for more');
    
end


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

pairs = isolateTrials(data, fullComparison);

%% Remove 'pairs' cells that are duplicates

pairs = dupeCellDelete(pairs);


%% Save data


%get filename
saveMatFile = input('choose a name for the full data file\n(optional, default used if left blank): ', 's');

%if no name is given, make one up using the current clock time and save
%a .mat file anyway
if length(saveMatFile)<1
    saveMatFile = sprintf('trimmedData_%s', datestr(now,'HHMM'));
    fprintf('\n(I am going to save your data anyway...\nit is called %s.mat)\n\n', saveMatFile);
end

%<CSV & DAT>
%save to a csv file
writeTrimmedData(pairs, saveMatFile);

%<WORKSPACE>
%place results in variables that will be available after the
%function has run.
varData = 'trimmedTrialsPerCell';

%assign variable into the specified name variable
assignin('base', varData, pairs);

%<MAT>
%save the workspace variables to a .MAT file 
save(saveMatFile);

%% End

fprintf('Finished. %s.dat is saved in this directory and ready for importing\ninto the MultiMatch function\n\n', saveMatFile);
