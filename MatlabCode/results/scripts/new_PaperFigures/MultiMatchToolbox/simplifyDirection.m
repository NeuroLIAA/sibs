function spGlobal = simplifyDirection(sp,T,Tdur)
%--------------------------------------------------------------------------
% Merges a series of short saccades
% 
% Input:    sp - scanpath vectors
%           T - threshold separating local from global scanpath
% Output:   sp_out - simplified scanpath
%
% Suggested modifications: Change 'and' to 'or'?
%--------------------------------------------------------------------------

%predefine variables
spGlobal.fixation = [];
spGlobal.saccade = [];

%--------------------------------------------------------------------------
% Calculate new vectors from consecutive local saccades
%--------------------------------------------------------------------------

if length(sp.saccade.x) <=1 %if there is only one or no starting x positions of a saccade recorded
    spGlobal = sp; %it cannot be simplified, so just use what's already assigned to sp
    return %end function
end

%set counters
i = 1;
j = 1;

% spGlobal.fixation.mu = sp.fixation.mu;
durMem = 0; %Presumably written in case the keepSaccade would track duration separately?
%in any case, as far as I can see it will always = 0.

% Can't have zero saccades...
while i <= length(sp.saccade.x)
    
    % Angular difference between two consecutive saccades
    if i < length(sp.saccade.x)
        v1=[sp.saccade.lenx(i),sp.saccade.leny(i)]; %(x,y) starting pos vector
        v2=[sp.saccade.lenx(i+1),sp.saccade.leny(i+1)]; %(x,y) end pos vector
        angle = acosd(dot(v1,v2)/(norm(v1)*norm(v2))); %calculate angle
    else
        angle = inf; %if this is the last entry, no angle.
    end

    % If the angle between two consecutive vectors is smaller than the threshold T, merge 
    % saccades i and i+1 (i.e. don't consider them to be separate saccades)
    if angle < T && i < length(sp.saccade.x)
        
        % Do not merge saccades if the intermediate fixation durations are
        % long (e.g. two saccades downwards would be merged even though
        % they would clearly be broken up by fixations)
        if sp.fixation.dur(i+1) >= Tdur
            [sp,spGlobal,i,durMem] = keepSaccade(sp,spGlobal,i,j,durMem);
            j = j+1;
            continue, %skip the remainder of the code inside this 'if' statement
        end
            
        % calculate sum of local vectors.
        v_x = sp.saccade.lenx(i) + sp.saccade.lenx(i+1);
        v_y = sp.saccade.leny(i) + sp.saccade.leny(i+1);
        [theta,len] = cart2pol(v_x,v_y);

        % ... save them a new global vectors
        spGlobal.saccade.x(j) = sp.saccade.x(i);
        spGlobal.saccade.y(j) = sp.saccade.y(i);        
        spGlobal.saccade.lenx(j) = v_x;
        spGlobal.saccade.leny(j) = v_y;
        spGlobal.saccade.len(j) = len;   
        spGlobal.saccade.theta(j) = theta;    
    
        %... and sum up all the fixation durations
        spGlobal.fixation.dur(j) = sp.fixation.dur(i);%+sp.fixation.dur(i+1)/2+durMem;
        durMem = 0;%sp.fixation.dur(i+1)/2;
        i = i+2;
        
    % If the last saccade in a scanpath has an amplitude < T
    elseif angle < T && i == length(sp.saccade.x)
        
        % Merge saccades only if the fixation durations are short        
        if sp.fixation.dur(i+1) >= Tdur
            [sp,spGlobal,i,durMem] = keepSaccade(sp,spGlobal,i,j,durMem);
            j = j+1;
            continue,
        end
        
        % calculate sum of local vectors.
        v_x = sp.saccade.lenx(i-1) + sp.saccade.lenx(i);
        v_y = sp.saccade.leny(i-1) + sp.saccade.leny(i);
        [theta,len] = cart2pol(v_x,v_y);

        % ... save them a new global vectors
        spGlobal.saccade.lenx(j-1) = v_x;
        spGlobal.saccade.leny(j-1) = v_y;
        spGlobal.saccade.len(j-1) = len;   
        spGlobal.saccade.theta(j-1) = theta;    
        spGlobal.fixation.dur(j-1) = spGlobal.fixation.dur(j-1)+sp.fixation.dur(i)/2;
%         spGlobal.fixation.dur = spGlobal.fixation.dur(1:end-1);
        durMem = 0;%sp.fixation.dur(i)/2;
        
        j = j-1;
        i = i+1;
        
    % Else just copy the saccade
    else
        [sp,spGlobal,i,durMem] = keepSaccade(sp,spGlobal,i,j,durMem);
    end
%     i,j
    j = j+1;
end

spGlobal.fixation.dur(j) = sp.fixation.dur(i)+durMem;


%---------------------------------------------------------------------
function [sp,spGlobal,i,durMem] = keepSaccade(sp,spGlobal,i,j,durMem)
%---------------------------------------------------------------------

spGlobal.saccade.x(j) = sp.saccade.x(i);
spGlobal.saccade.y(j) = sp.saccade.y(i);
spGlobal.saccade.lenx(j) = sp.saccade.lenx(i);
spGlobal.saccade.leny(j) = sp.saccade.leny(i);
spGlobal.saccade.len(j) = sp.saccade.len(i);
spGlobal.saccade.theta(j) = sp.saccade.theta(i);
spGlobal.fixation.dur(j) = sp.fixation.dur(i)+durMem;
durMem = 0;

i = i+1;




