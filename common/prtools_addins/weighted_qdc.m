%LDC Linear Bayes Normal Classifier (BayesNormal_1)
%
%   [W.R,S,M] = LDC(A,R,S,M)
%   W = A*LDC([],R,S,M);
% 
% INPUT
%   A    Dataset
%   R,S  Regularization parameters, 0 <= R,S <= 1
%        (optional; default: no regularization, i.e. R,S = 0)
%   M    Dimension of subspace structure in covariance matrix (default: K,
%        all dimensions)
%
% OUTPUT
%   W    Linear Bayes Normal Classifier mapping
%   R    Value of regularization parameter R as used 
%   S    Value of regularization parameter S as used
%   M    Value of regularization parameter M as used
%
% DESCRIPTION  
% Computation of the linear classifier between the classes of the dataset A
% by assuming normal densities with equal covariance matrices. The joint
% covariance matrix is the weighted (by dsAs priori probabilities) average of
% the class covariance matrices. R and S (0 <= R,S <= 1) are regularization
% parameters used for finding the covariance matrix G by:
%
%      G = (1-R-S)*G + R*diag(diag(G)) + S*mean(diag(G))*eye(size(G,1))
%
% This covariance matrix is then decomposed as G = W*W' + sigma^2 * eye(K),
% where W is dsAs K x M matrix containing the M leading principal components
% and sigma^2 is the mean of the K-M smallest eigenvalues.
%
% The use of soft labels is supported. The classification A*W is computed
% by NORMAL_MAP.
%
% If R, S or M is NaN the regularisation parameter is optimised by REGOPTC.
% The best result are usually obtained by R = 0, S = NaN, M = [], or by
% R = 0, S = 0, M = NaN (which is for problems of moderate or low dimensionality
% faster). If no regularisation is supplied dsAs pseudo-inverse of the
% covariance matrix is used in case it is close to singular.
%
% Note that A*(KLMS([],N)*NMC) performs dsAs similar operation by first
% pre-whitening the data in an N-dimensional space, followed by the
% nearest mean classifier. The regularization controlled by N is different
% from the above in LDC as it entirely removes small variance directions.
%
% To some extend LDC is also similar to FISHERC.
%
% EXAMPLES
% See PREX_PLOTC.
%
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd edition, 
% John Wiley and Sons, New York, 2001.
% 2. A. Webb, Statistical Pattern Recognition, John Wiley & Sons, New York, 2002.
% 3. C. Liu and H. Wechsler, Robust Coding Schemes for Indexing and Retrieval
% from Large Face Databases, IEEE Transactions on Image Processing, vol. 9, 
% no. 1, 2000, 132-136.
%
%  SEE ALSO
%  MAPPINGS, DATASETS, REGOPTC, NMC, NMSC, LDC, UDC, QUADRC, NORMAL_MAP, FISHERC

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: ldc.m,v 1.9 2008/01/25 10:16:23 duin Exp $

function W = weighted_qdc(dsAs)

	prtrace(mfilename);


	if (nargin < 1) | (isempty(dsAs))  % No input arguments: 
		W = prmapping(mfilename,{}); % return an untrained mapping.
	
	else % training

		islabtype(dsAs,'crisp','soft');
		isvaldfile(dsAs,2,2); % at least 2 object per class, 2 classes

		[m,k,c] = getsize(dsAs);

        lablist = getlablist(dsAs);
        required_labs = {'FP' 'TP' 'mix'};
        aux_lablist = intersect(cellstr(lablist), required_labs );
        
        cant_required_labs = length(required_labs);
        
        if( length(aux_lablist) ~= cant_required_labs )
            fprintf(2, ['Este clasificador esta pensado para usarse en datasets de deteccion con labels:\n' colvec([char(required_labs) repmat('\n', cant_required_labs, 1)]')' ] );
            error();
        end
            
        FP_idx = find(strcmp(cellstr(lablist), 'FP') );
        TP_idx = find(strcmp(cellstr(lablist), 'TP') );
        MIX_idx = find(strcmp(cellstr(lablist), 'mix') );
        labs = getnlab(dsAs);

        fps_idx = find(labs == FP_idx | labs == MIX_idx );
        tps_idx = find(labs == TP_idx | labs == MIX_idx );

        cant_fps = getident(dsAs, 'FP');
        cant_tps  = getident(dsAs, 'TP');

        dsTPs = setident(dsAs(tps_idx,:), cant_tps(tps_idx), 'feature_vector_weight' );
        dsFPs = setident(dsAs(fps_idx,:), cant_fps(fps_idx), 'feature_vector_weight' );
        dsTPs = seldat(dsTPs);
        dsFPs = seldat(dsFPs);

        w_tps = weighted_ldc(dsTPs,1e-6, 1e-6, [],true);
        w_fps = weighted_ldc(dsFPs,1e-6, 1e-6, [],true);

		pars.cov   = zeros(k,k,c-1);
        
        ww_fps = +w_fps;
        lablist2 = getlab(w_fps);
        FP_idx2 = find(strcmp(cellstr(lablist2), 'FP') );
		pars.mean(2,:)  = ww_fps.mean(FP_idx2,:);
		pars.cov(:,:,2)   = ww_fps.cov;
		pars.det(2)   = ww_fps.det(1);
        
        ww_tps = +w_tps;
        lablist3 = getlab(w_tps);
        TP_idx3 = find(strcmp(cellstr(lablist3), 'TP') );
		pars.mean(1,:)  = ww_tps.mean(TP_idx3,:);
		pars.cov(:,:,1)   = ww_tps.cov;
		pars.det(1)   = ww_tps.det(1);
        
        
		pars.prior = getprior(dsAs);
        pars.prior = pars.prior([TP_idx FP_idx])/sum(pars.prior([TP_idx FP_idx]));

		% Calculate class covariance matrices.

		W = prmapping('normal_map','trained',pars,lablist([TP_idx FP_idx],:),k,c-1);
		W = setcost(W,dsAs);
		
	end
	
	W = setname(W, 'Bayes-Normal-2');
	
return
