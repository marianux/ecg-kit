function [autovec autoval]= autovec_calculation_robust(wtECGslice)

    mean_wtECGslice = mean(wtECGslice);    
    result = mcdcov( bsxfun( @minus, wtECGslice, mean_wtECGslice ), 'plots', 0);
    wtECGslice_cov = result.cov;
    [autovec autoval] = eig(wtECGslice_cov); 
    autoval = diag(autoval);
    [~, autoval_idx] = sort(autoval, 'descend');
    autovec = autovec(:,autoval_idx);
    autoval = autoval(autoval_idx);
