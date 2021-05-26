function sp1 = generateStructureArrayScanpath(data)

%obtain number of rows (LENgth)
len = size(data,1);

%for every row
for i = 1:len
    
    %put the x, y and duration positions of each row into a structure
    %variable
    sp1.fixation.x(i) = data(i,1);         
    sp1.fixation.y(i) = data(i,2);                         
    sp1.fixation.dur(i) = data(i,3);   
    
    %if i is less than the number of rows (i.e. perform all lines inside
    %this if statement unless this is the last row of the data)
    if i < len
        sp1.saccade.x(i) = data(i,1); %starting x position of a saccade   
        sp1.saccade.y(i) = data(i,2); %starting y position of a saccade   
    end
    
    %run lines inside this if statement as long as this is not the first
    %fixation
    if i > 1
        
        %saccade x-distance travelled
        sp1.saccade.lenx(i-1) = (sp1.fixation.x(i) - sp1.saccade.x(i-1));
        
        %saccade y-distance travelled
        sp1.saccade.leny(i-1) = (sp1.fixation.y(i) - sp1.saccade.y(i-1));
        
        %obtain angle (theta) and radius (len) of the saccade,
        %respectively.
        [sp1.saccade.theta(i-1),sp1.saccade.len(i-1)] = cart2pol(sp1.saccade.lenx(i-1),sp1.saccade.leny(i-1));
    end

    
end

