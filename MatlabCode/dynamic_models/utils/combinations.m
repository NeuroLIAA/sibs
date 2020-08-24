function result=combinations(elements)
    combinations = cell(1, numel(elements)); %set up the varargout result
    [combinations{:}] = ndgrid(elements{:});
    combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); %there may be a better way to do this
    result = [combinations{:}]; % NumberOfCombinations by N matrix. Each row is unique.
end