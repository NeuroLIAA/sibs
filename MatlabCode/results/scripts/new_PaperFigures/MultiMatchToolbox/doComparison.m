function rv =doComparison(data1,data2, sz, suppressOutput, saveName)
%rv = doComparison(fixationData1,fixationData2[, sz, saveName])
%
% Created by Tom Foulsham (foulsham@essex.ac.uk), modified by Chris Street
% (c.street@ucl.ac.uk)
%
% [This version of the script has been edited by Chris Street
% (c.street@ucl.ac.uk/cstreet1986@gmail.com) on 17/09/12 to improve
% commenting, layout style and script performance]
%
%Pass in two n x 3 fixation arrays and do a MultiMatch comparison. Each
%array should be organised so that the first column is x position, the
%second y position, and the third duration of a given fixation, with each
%row representing a fixation for a given trial. The sz
%input is the size of the screen the data were collected on, in the form
%[x y]. This is an optional input, but will be requested by the script 
%if left blank. Finally, a name (which must be inputted as a string) can be entered
%to save all the variables into a .mat file. This is optional.
%
%returns rv, a 5 item row vector with similarity values for each
%dimension, which are: (1) euclidean difference, (2) direction difference,
%(3) length difference, (4) position difference, and (5) duration difference
%between a given fixation from one matrix (data1) and the matching fixation 
%in the other matrix (data2)
%
%rv is NaN when ALL the fixations were off screen during stim presentation.


% sz = [1024 768];%adjust parameters for JOV experiment


%obtain screen size
if nargin==2
    sz = input('screen size (e.g. [800 600]): ');
end

%check screen size makes sense
if length(sz) ~= 2 %if not equal to 2
    error('[x y] screen size must be specified');
end
    

%CHANGE GLOBAL THRESH ACCORDINGLY?
%set thresholds for detecting saccades (globalThreshold), determining
%direction of saccade and (the duration of a saccade?)
globalThreshold = 0.1*sqrt(sz(1)^2+sz(2)^2); % Currently, screen diagonal/10.
directionThreshold = 45;
durationThreshold = inf; % (what does this do?)

%return NaN  and error if scanpaths do not have 3 columns or at least 2 rows
if size(data1,2) ~= 3 || size(data2,2) ~=3 || size(data1,1) < 2  ...
        || size(data2,1) < 2 
    rv=NaN;
%    return
    error(sprintf('Each inputted array must contain two fixations and three columns\nthat specify details of the fixation (x,y,duration)\n'));

else %if data are acceptable for use
    
%     bx1=data1(:,1)>0 & data1(:,1) <1024;
%     bx2=data2(:,1)>0 & data2(:,1) <1024;
%     by1=data1(:,2)>0 & data1(:,2) <768;
%     by2=data2(:,2)>0 & data2(:,2) <768;
    
    %get data from the two arrays that are smaller than the confines of 
    %the screen (greater than 0, less than screen size)
    bx1=data1(:,1)>0 & data1(:,1) <sz(1); 
    bx2=data2(:,1)>0 & data2(:,1) <sz(1);
    by1=data1(:,2)>0 & data1(:,2) <sz(2);
    by2=data2(:,2)>0 & data2(:,2) <sz(2);
    

    %return NaN if any fixations go out of bounds
    %TOO CONSERVATIVE? (but none of them will because the above data
    %sifting removed any values that went out of bounds)
    %i.e. only perform enclosed lines if all the x,y positions are nonzero,
    %as they should be as a result of the previous lines that filtered the
    %code
    if all(bx1) && all(bx2) && all(by1) && all(by2)
    
        %transform into scanpath structure
        sp1 = generateStructureArrayScanpath(data1);
        sp2 = generateStructureArrayScanpath(data2);

        %--------------------------------------------------------------------------
        % Compare scanpaths 
        %--------------------------------------------------------------------------

        % ... using the proposed method
        [sp1,sp2,rv,path,M_assignment,ki,kj] = mainProposed(sp1,sp2,...
            globalThreshold,directionThreshold,durationThreshold,sz);

        if suppressOutput == 0 %if no save name given, best display the data
            %print results to the prompt
            disp(['Vector similarity = ',num2str(rv(1))])
            disp(['Direction similarity = ',num2str(rv(2))])
            disp(['Length similarity = ',num2str(rv(3))])
            disp(['Position similarity = ',num2str(rv(4))])
            disp(['Duration similarity = ',num2str(rv(5))])
        end
        
        %place results in variables that will be available after the
        %function has run.
        varVec = 'vectorSim';
        varDir = 'directionSim';
        varLength = 'lengthSim';
        varPos = 'positionSim';
        varDur = 'durationSim';
        scanpath1 = 'scanpath1';
        scanpath2 = 'scanpath2';
        
        %assign variables into the specified name variables
        assignin('base', varVec, rv(1));
        assignin('base', varDir, rv(2));
        assignin('base', varLength, rv(3));
        assignin('base', varPos, rv(4));
        assignin('base', varDur, rv(5));
        assignin('base', scanpath1, sp1);
        assignin('base', scanpath2, sp2);
        
        if nargin == 5 %if a saveName is specified
           
            %ensure input is a valid variable name to be used - if not,
            %generate one
            varName = genvarname(saveName);
            saveName = varName; %duplicate varName
            
            %check saveName ends with ".mat" - if it doesnt, make it end with
            %that extension
            extn = saveName(length(saveName)-3:length(saveName));
            if strcmp(extn, '.mat') == 0 %if it doesn't end with .mat
                saveName = strcat(saveName, '.mat');
            end
            
            %assign variable into the specified name variable
            save(varName);
            
            %****move data to the /mat_files subfolder*****
            %create subdirectory
            if ~exist('mat_files', 'dir')
                mkdir('mat_files');
            end
            
            %get destination as string
            curdir = cd;
            if ismac == 1
                dest = strcat(curdir, '/mat_files'); %get string for wanted directory
            else
                dest = strcat(curdir, '\mat_files'); %get string for wanted directory
            end
            
            %move to destination
            movefile(saveName, dest, 'f');
            
        end
    
    else %if any fixation position is 0 (i.e. on boundary of the screen)
        
        %output no value and error.
        rv=NaN;
        
        %show reason for error if output hasn't been suppressed
        if suppressOutput == 0
            s1 = ('a fixation was either on or outwith the screen boundaries. This code\nshould filter out any of these cases before');
            s2 = (' running comparisons,\nbut has failed to do so. This is likely to be a coding error, but removing any\n');
            s3 = ('(x,y) positions that would be at or over the screen boundary may work\n\n');
            %concatenate above strings
            errtext = strcat(s1, s2, s3);
    %         %display as error
    %         error(sprintf(errtext));

            %display problem to screen
            fprintf(errtext);
        end
        
    end
end




