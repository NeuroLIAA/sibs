function result = fun_heatmap(im_file, gaze_x, gaze_y, settings, name)
%genera heatmap. las coordenadas del gaze tienen que estar centradas en la
%figura
%
% requiere fun_hist2d.m
%
%ultimo update: Jul 29 2016 (GS)


%% importa imagen
im          = imread(im_file);
if settings.grayscale
    if size(size(im))==3
        im = rgb2gray(im);
    else
        im = repmat(im, [1 1 3]);
    end
end

[ny, nx, ~] = size(im);

%% construyo histograma suavizado
nhist = fun_hist2d(gaze_x, gaze_y, 1:nx, 1:ny);

gfilt = fspecial('gaussian', settings.filterSize, settings.filterSD);
nhist = imfilter(nhist,gfilt);
nhist = nhist /max(max(nhist));
nhist(nhist<0)=0;
nhist = nhist.^.4;

%%
result = settings.alpha * nhist;
if nargin == 5
    imwrite(settings.alpha * nhist, name);
    return;
end

if nargout == 0
    imshow(im);
    hold on
    h = image(settings.max_nhist * nhist);
    colormap(jet);
    set(h, 'AlphaData', settings.alpha * nhist)
end
end