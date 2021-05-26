function models = fun_define_models()
    j=1;
        models(j).prior     = 'deepgaze';
        models(j).searcher  = 'correlation';
        models(j).params    = 'a_3_b_4_tam_celda_32';
        models(j).cols      = 'r';
        models(j).name      = 'dg:cibs';

    j=2;
        models(j).prior     = 'deepgaze';
        models(j).searcher  = 'geisler';
        models(j).params    = 'a_3_b_4_tam_celda_32';
        models(j).cols      = 'g';
        models(j).name      = 'dg:ibs';

    j=3;
        models(j).prior     = 'deepgaze';
        models(j).searcher  = 'greedy';
        models(j).params    = 'a_3_b_4_tam_celda_32';
        models(j).cols      = 'b';
        models(j).name      = 'dg:greedy';

    j=4;
        models(j).prior     = 'center';
        models(j).searcher  = 'correlation';
        models(j).params    = 'a_3_b_4_tam_celda_32';
        models(j).cols      = 'y';
        models(j).name      = 'ctr:cibs';

    j=5;
        models(j).prior     = 'noisy';
        models(j).searcher  = 'correlation';
        models(j).params    = 'a_3_b_4_tam_celda_32';
        models(j).cols      = 'y';
        models(j).name      = 'noisy:cibs';

    j=6;
        models(j).prior     = 'flat';
        models(j).searcher  = 'correlation';
        models(j).params    = 'a_3_b_4_tam_celda_32';
        models(j).cols      = 'y';
        models(j).name      = 'flat:cibs';

    j=7;
        models(j).prior     = 'mlnet';
        models(j).searcher  = 'correlation';
        models(j).params    = 'a_3_b_4_tam_celda_32';
        models(j).cols      = 'r';
        models(j).name      = 'mlnet:cibs';

    j=8;
        models(j).prior     = 'sam-resnet';
        models(j).searcher  = 'correlation';
        models(j).params    = 'a_3_b_4_tam_celda_32';
        models(j).cols      = 'r';
        models(j).name      = 'sam-resnet:cibs';

    j=9;
        models(j).prior     = 'sam-vgg';
        models(j).searcher  = 'correlation';
        models(j).params    = 'a_3_b_4_tam_celda_32';
        models(j).cols      = 'r';
        models(j).name      = 'sam-vgg:cibs';
end