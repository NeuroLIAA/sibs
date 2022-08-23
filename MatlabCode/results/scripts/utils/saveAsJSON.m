function saveAsJSON(data, jsonFileName)
% saves the values in the structure 'data' to a file in JSON format.
% Based on the work of Lior Kirsch at: https://uk.mathworks.com/matlabcentr
% al/fileexchange/50965-structure-to-json
% 
% Modification by Arthur Y.C. Liu 24/06/2020
% 
% Example:
%     data.name = 'chair';
%     data.color = 'pink';
%     data.eye = eye(3);
%     data.metrics.imageSize = [1920; 1080];
%     data.metrics.height = 0.3;
%     data.metrics.width = 1.3;
%     saveJSONfile(data, 'out.json');
% 
% Output 'out.json':
% {
% 	"name" : "chair",
% 	"color" : "pink",
% 	"eye" : 
% 	[
% 		1,
% 		0,
% 		0,
% 		0,
% 		1,
% 		0,
% 		0,
% 		0,
% 		1
% 	],
% 	"metrics" : 
% 	{
% 		"imageSize" : 
% 		[
% 			1920,
% 			1080
% 		],
% 		"height" : 0.3,
% 		"width" : 1.3
% 	}
% }
fid = fopen(jsonFileName,'w');
if isobject(data)
    data = toStruct(data);
end
writeElement(fid, data,'', true);
fprintf(fid,'\n');
fclose(fid);
end
function writeElement(fid, data, tabs, isFirstLine)
namesOfFields = fieldnames(data);
numFields = length(namesOfFields);
tabs = sprintf('%s\t', tabs);
if nargin == 4
    if isFirstLine
        fprintf(fid,'%s{\n%s', tabs(1:end-1), tabs);
    end
else
    fprintf(fid,'\n%s{\n%s', tabs(1:end-1), tabs);
end
for i = 1:numFields
    currentField = namesOfFields{i};
    currentElementValue = data.(currentField);
    writeSingleElement(fid, currentField, currentElementValue, tabs);
    if i == numFields
        fprintf(fid,'\n%s}',  tabs(1:end-1)); 
    else
        fprintf(fid,',\n%s', tabs);
    end
end
end
function writeSingleElement(fid, currentField, currentElementValue, tabs)
% if this is an array/matrix and not a string then iterate on every
% element, if this is a single element write it
if ~isstruct(currentElementValue) &&...
        length(currentElementValue) > 1 && ~ischar(currentElementValue)  
    fprintf(fid,'"%s" : \n%s[\n%s\t',currentField, tabs, tabs);
    tabs = sprintf('%s\t', tabs);
    valLength = length(currentElementValue(:));
    for m = 1:valLength
        fprintf(fid,'%g' , currentElementValue(m));
        if m == valLength
            fprintf(fid,'\n%s]', tabs(1:end-1));
        else
            fprintf(fid,',\n%s', tabs);
        end
    end    
elseif isstruct(currentElementValue) &&...
        length(currentElementValue) > 1 && ~ischar(currentElementValue)
    fprintf(fid,'"%s" : \n%s[',currentField, tabs);
    tabs = sprintf('%s\t', tabs);
    valLength = length(currentElementValue(:));
    for m = 1:valLength
        writeElement(fid, currentElementValue(m), tabs);
        if m == valLength
            fprintf(fid,'\n%s]', tabs(1:end-1));
        else
            fprintf(fid,',');
        end
    end
elseif isstruct(currentElementValue)
    fprintf(fid,'"%s" : ',currentField);
    writeElement(fid, currentElementValue, tabs);
elseif isnumeric(currentElementValue) || islogical(currentElementValue)
    fprintf(fid,'"%s" : %g' , currentField, currentElementValue);
elseif isempty(currentElementValue)
    fprintf(fid,'"%s" : "null"' , currentField);
elseif isobject(currentElementValue)
    currentElementValue = toStruct(currentElementValue);
    fprintf(fid,'"%s" : ',currentField);
    writeElement(fid, currentElementValue, tabs);
else %ischar or something else ...
    fprintf(fid,'"%s" : "%s"' , currentField, currentElementValue);
end
end