myFiles = dir('profiles/*.mat'); %gets all wav files in struct
for k = 1:length(myFiles)
    load(fullfile(myFiles(k).folder, myFiles(k).name))
    display(myFiles(k).name)
    display(max([p.FunctionTable(:).TotalTime]))
end
