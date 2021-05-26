function [sumdata sumdataString participantMerge] = readFromHeadedFile(...
            saveMatFile)
%will read data from a .csv file that contains a header row. Dismisses that
%row and leaves you with the data in a form that can be used by the
%createParticipantSummaries() function in MultiMatch. Heavily reliant on work posted in 
%http://stackoverflow.com/questions/2588021/skip-reading-headers-in-matlab 
        

%% sumdata and sumdataString
%get filename
fileName = strcat(saveMatFile, '.csv');
fid = fopen(fileName); %open filename

% read all text (15 columns)
strData = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','Delimiter',','); 

%close the file again
fclose(fid);

sumdataStringHeader = strData{5};
sumdataString = sumdataStringHeader(2:end); %remove header <<< SUMDATASTRING

%concatenate because textscan returns a column of cells for each column in the data
strData = cat(2,strData{:}); 

%# convert columns 1:4 to double
doubleDataA = str2double(strData(:,1:4));

%# convert columns 7:15 to double
doubleDataB = str2double(strData(:,7:end));

%concatenate doubleData
doubleData = [doubleDataA doubleDataB];

%# find header rows. headerRows is a logical array
headerRowsL = all(isnan(doubleData),2);

%remove headers, results in sumdata  <<< SUMDATA
sumdataTemp = doubleData(~headerRowsL,:);
%add two columns of zeros in the middle
sumdata = [sumdataTemp(:,1:4) zeros(size(sumdataTemp,1),2) sumdataTemp(:,5:end)];


%% participantMerge

%get a concatenation of the two participants
participantMerge = textread('pMerge.csv', '%s'); % <<<PARTICIPANTMERGE
