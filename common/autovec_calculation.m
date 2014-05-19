function [autovec autoval]= autovec_calculation(wtECGslice)

    mean_wtECGslice = mean(wtECGslice);    
    wtECGslice_cov = cov( bsxfun( @minus, wtECGslice, mean_wtECGslice ));
    [autovec autoval] = eig(wtECGslice_cov); 
    autoval = diag(autoval);
    [~, autoval_idx] = sort(autoval, 'descend');
    autovec = autovec(:,autoval_idx);
    autoval = autoval(autoval_idx);
