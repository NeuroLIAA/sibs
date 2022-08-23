clear all
close all
clc
a = dir('../images/*.jpg'); imgnames = {a.name};
img_list = [];
for i = 116:134%length(imgnames) 
    img = imread(['../images/' imgnames{i}]); %figure(1); clf; imagesc(img)
    tmp = dir(['../templates/' imgnames{i}(1:end-4) '*']);
    template = imread(['../templates/' tmp.name]); 
    target_size = size(template);

    %ssimval = calculate_ssim(img,template);
    %if ~exist('/ssim_precalc'); mkdir('/ssim_precalc'); end
    %save(['./ssim_precalc/ssimval_pad_' imgnames{i}(1:end-4)],'ssimval');
end