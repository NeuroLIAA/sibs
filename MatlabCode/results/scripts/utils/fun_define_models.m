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
        case 'priors-ssim' % deepgaze searchers + saliency
            disp('Loading priors-sIBS...')
            j = 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'structuralsim';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.8 0.6 0.2];
                models(j).name          = 'sIBS+DGII';
            j = j + 1;
                models(j).prior         = 'center'; 
                models(j).searcher      = 'structuralsim';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.5 0.8 0.2];
                models(j).name          = 'sIBS+Center';
            j = j + 1;
                models(j).prior         = 'noisy'; 
                models(j).searcher      = 'structuralsim';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.6 0.1 0.8];
                models(j).name          = 'sIBS+Noisy';
            j = j + 1;
                models(j).prior         = 'flat'; 
                models(j).searcher      = 'structuralsim';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.1 0.6 0.5];
                models(j).name          = 'sIBS+Flat';
        case 'priors-correlation'  % cIBS priors
            disp('Loading priors-cIBS...')
            j = 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'correlation';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'r';
                models(j).name          = 'cIBS+DGII';
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
        case 'priors-ibs' %
            disp('Loading priors-IBS...')
            j = 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'b';
                models(j).name          = 'IBS+DGII';
            j = j + 1;
                models(j).prior         = 'center'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.9 0.3 0.4];
                models(j).name          = 'IBS+Center';
            j = j + 1;
                models(j).prior         = 'flat'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.6 0.5 0.1];
                models(j).name          = 'IBS+Flat';
            j = j + 1;
                models(j).prior         = 'noisy'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.2 0.8 0.4];
                models(j).name          = 'IBS+Noisy';
        case '5searchers'%searchers-deepgaze-ssim' % deepgaze searchers + saliency
            j = 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'correlation';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'r';
                models(j).name          = 'cIBS+DGII';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'structuralsim';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.8 0.6 0.2];
                models(j).name          = 'sIBS+DGII';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'b';
                models(j).name          = 'IBS+DGII';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'greedy';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'g';
                models(j).name          = 'Greedy+DGII';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'saliency-based';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'm';
                models(j).name          = 'SaliencyBased';
        case 'all'
            j = 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'correlation';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'r';
                models(j).name          = 'cIBS+DGII';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'b';
                models(j).name          = 'IBS+DGII';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'greedy';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'g';
                models(j).name          = 'Greedy+DGII';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'saliency-based';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'm';
                models(j).name          = 'SaliencyB.';
            j = j + 1;
                models(j).prior         = 'flat'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'r';
                models(j).name          = 'IBS+Flat';
            j = j + 1;
                models(j).prior         = 'center'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'b';
                models(j).name          = 'IBS+Center';
            j = j + 1;
                models(j).prior         = 'noisy'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'b';
                models(j).name          = 'IBS+Noisy';
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
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'structuralsim';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.8 0.6 0.2];
                models(j).name          = 'sIBS+DGII';
            j = j + 1;
                models(j).prior         = 'center'; 
                models(j).searcher      = 'structuralsim';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.5 0.8 0.2];
                models(j).name          = 'sIBS+Center';
            j = j + 1;
                models(j).prior         = 'noisy'; 
                models(j).searcher      = 'structuralsim';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.6 0.1 0.8];
                models(j).name          = 'sIBS+Noisy';
            j = j + 1;
                models(j).prior         = 'flat'; 
                models(j).searcher      = 'structuralsim';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = [0.1 0.6 0.5];
                models(j).name          = 'sIBS+Flat';
        case 'searchers-deepgaze' % deepgaze searchers + saliency
            j = 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'correlation';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'r';
                models(j).name          = 'cIBS+DGII';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'geisler';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'b';
                models(j).name          = 'IBS+DGII';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'greedy';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'g';
                models(j).name          = 'Greedy+DGII';
            j = j + 1;
                models(j).prior         = 'deepgaze'; 
                models(j).searcher      = 'saliency-based';
                models(j).params        = 'a_3_b_4_tam_celda_32';
                models(j).cols          = 'm';
                models(j).name          = 'SaliencyBased';
    end
end
