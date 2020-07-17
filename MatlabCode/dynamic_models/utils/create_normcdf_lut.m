function r = create_normcdf_lut(e)
    x = norminv(e:e:1-e);
    y = normcdf(x);
    tab.x = x;
    tab.y = y;
    r = tab;
   