%% (Internal) Calculate eigenvalues and vectors
%   
%   [autovec autoval]= autovec_calculation(wtECGslice)
% 
% Arguments:
% 
%      + wtECGslice: the multivariate signal
% 
% Output:
% 
%      + autovec: eigenvectors
% 
%      + autoval: eigenvalues
% 
% Example:
% 
%     % 1,2 PCA at scale 4
%     rotation_matrix = autovec_calculation(wt_mean_hb_pos_cell{scale_idx(4)});
%     rotation_matrix = rotation_matrix(:,1:2);
% 
% 
% See also BaselineWanderRemovalMedian
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function [autovec autoval]= autovec_calculation(wtECGslice)

    mean_wtECGslice = mean(wtECGslice);    
    wtECGslice_cov = cov( bsxfun( @minus, wtECGslice, mean_wtECGslice ));
    [autovec autoval] = eig(wtECGslice_cov); 
    autoval = diag(autoval);
    [~, autoval_idx] = sort(autoval, 'descend');
    autovec = autovec(:,autoval_idx);
    autoval = autoval(autoval_idx);
