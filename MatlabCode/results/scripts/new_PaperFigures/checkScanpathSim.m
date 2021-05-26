function ssim = checkScanpathSim(index, scan, grid_size, infoc)
    if length(infoc)==1
        index = 1;
    end
    ssim = scanpathDistance([infoc(index).fixations_matrix_reduced(:,2), infoc(index).fixations_matrix_reduced(:,1)], ...
                scan, grid_size);
    
end