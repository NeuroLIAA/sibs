function auc = fun_auc_with_errors(meanSalientPercent,meanTPR,errTPR)
    auc.mean    = trapz(meanSalientPercent, meanTPR);               
    auc.err(1)	= trapz(meanSalientPercent, meanTPR-errTPR);             
    auc.err(2)	= trapz(meanSalientPercent, meanTPR+errTPR);             
end