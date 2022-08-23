clear all
close all
clc

%%

a = dir('../images/*.jpg'); imgnames = {a.name};
b = dir('./ssim_precalc/*.mat'); ssimnames = {b.name}';
prefix = 'ssimval_pad_';
suffix = '.mat';
ssimnames = cellfun(@(x) strrep(x,prefix,''),ssimnames,'UniformOutput', false);
ssimnames = cellfun(@(x) strrep(x,suffix,'.jpg'),ssimnames,'UniformOutput', false);

%%
for i = 1:85
    if ~any(strcmp(ssimnames,imgnames{i})) 
   
        img = imread(['../images/' imgnames{i}]); %figure(1); clf; imagesc(img)
        tmp = dir(['../templates/' imgnames{i}(1:end-4) '*']);
        template = imread(['../templates/' tmp.name]); 
        target_size = size(template);
        ssimval = calculate_ssim(img,template);
        %if ~exist('/ssim_precalc'); mkdir('/ssim_precalc'); end
        save(['./ssim_precalc/ssimval_pad_' imgnames{i}(1:end-4)],'ssimval');
    else
        disp(imgnames{i})
    end
end