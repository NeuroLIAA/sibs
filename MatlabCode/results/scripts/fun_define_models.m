function models = fun_define_models(experiment)
% path = sprintf('../out_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
%                 models(ind_model).prior,...
%                 models(ind_model).searcher,...
%                 models(ind_model).params,...
%                 ind_img);
% Priors:
% 'noisy'
% 'flat'
% 'center' 
% 'icf' 
% 'deepgaze' 
% 'mlnet' 
% 'sam-resnet' 
% 'sam-vgg'
%
% Searchers:
% 'correlation'
% 'geisler'
% 'greedy'
%
% Params:
% 'a_3_b_4_tam_celda_32'

    switch experiment
        case 'searchers-deepgaze' % deepgaze searchers + saliency
            j = 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'correlation';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'r';
                models(j).name          = 'cIBS+DeepGaze';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'b';
                models(j).name          = 'IBS+DeepGaze';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'greedy';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'g';
                models(j).name          = 'Greedy+DeepGaze';
            j = j + 1;
                models(j).prior         = 'saliency-based'; 
                models(j).searcher      = 'correlation';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'm';
                models(j).name          = 'SaliencyBased';
        case 'priors-correlation'  % cIBS priors
            j = 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'correlation';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'r';
                models(j).name          = 'cIBS+DeepGaze';
            j = j + 1;
                models(j).prior         = 'center'; 
                models(j).searcher      = 'correlation';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0 0.5 1];
                models(j).name          = 'cIBS+Center';
            j = j + 1;
                models(j).prior         = 'flat'; 
                models(j).searcher      = 'correlation';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.5 0 1];
                models(j).name          = 'cIBS+Flat';
            j = j + 1;
                models(j).prior         = 'noisy'; 
                models(j).searcher      = 'correlation';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [1 0.6 0];
                models(j).name          = 'cIBS+Noisy';
        case 'ibs-flat-center' % Flat vs DG2
            j = 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'b';
                models(j).name          = 'IBS+DeepGaze';
            j = j + 1;
                models(j).prior         = 'flat'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'r';
                models(j).name          = 'IBS+Flat';
    end
end
