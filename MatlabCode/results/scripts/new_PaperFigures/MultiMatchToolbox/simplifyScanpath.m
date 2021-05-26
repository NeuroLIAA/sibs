function spGlobal = simplifyScanpath(spGlobal,T1,T2,Tdur)
%spGlobal = simplifyScanpath(spGlobal,T1,T2,Tdur);
%  inputs: sp = scanpath, T1 = global threshold, T2 = direction
%  threshold, Tdur = duration threshold. Outputs a simplified version of
%  the scanpath

l = inf; %predefine variable
% spGlobal = sp; %scanpath

% Loop until no further simplifications are made
% h = figure
while 1 %run continuously
    
    
    spGlobal = simplifyDirection(spGlobal,T2,Tdur);
%     spGlobal.fixation.dur
    spGlobal = simplifyLength(spGlobal,T1,Tdur);
%     spGlobal.fixation.dur

%     h = subplot(1,2,2)
%     plotScanpaths(spGlobal,spGlobal,'simplified sp',[1 1],[5 5])
%     drawnow
%     pause

    %if l (length of variableX at last run of this loop) = length of
    %variableX on this run of the loop (i.e. the size of the variable
    %hasn't increased), then stop running the loop
    
    if l == length(spGlobal.fixation.dur) %if the variable is infinitely long
        break %stop running
    end
    
    l = length(spGlobal.fixation.dur); %update length of variable
%     clf(h)
end

% close(h)



