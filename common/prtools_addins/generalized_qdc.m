%QDC Quadratic Bayes Normal Classifier (Bayes-Normal-2)
%
%   [W,R,S,M] = QDC(A,R,S,M)
%   W = A*QDC([],R,S)
%
% INPUT
%   A    Dataset
%   R,S	 Regularization parameters, 0 <= R,S <= 1 
%        (optional; default: no regularization, i.e. R,S = 0)
%   M    Dimension of subspace structure in covariance matrix (default: K,
%        all dimensions)
%
% OUTPUT
%   W    Quadratic Bayes Normal Classifier mapping
%   R    Value of regularization parameter R as used 
%   S    Value of regularization parameter S as used
%   M    Value of regularization parameter M as used
%
% DESCRIPTION
% Computation of the quadratic classifier between the classes of the dataset
% A assuming normal densities. R and S (0 <= R,S <= 1) are regularization
% parameters used for finding the covariance matrix by
% 
%   G = (1-R-S)*G + R*diag(diag(G)) + S*mean(diag(G))*eye(size(G,1))
%
% This covariance matrix is then decomposed as G = W*W' + sigma^2 * eye(K),
% where W is a K x M matrix containing the M leading principal components
% and sigma^2 is the mean of the K-M smallest eigenvalues.
%
% 
% 
% The use of soft labels is supported. The classification A*W is computed by
% NORMAL_MAP.
%
% If R, S or M is NaN the regularisation parameter is optimised by REGOPTC.
% The best result are usually obtained by R = 0, S = NaN, M = [], or by
% R = 0, S = 0, M = NaN (which is for problems of moderate or low dimensionality
% faster). If no regularisation is supplied a pseudo-inverse of the
% covariance matrix is used in case it is close to singular.
%
% EXAMPLES
% See PREX_MCPLOT, PREX_PLOTC.
%
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd
% edition, John Wiley and Sons, New York, 2001. 
% 2. A. Webb, Statistical Pattern Recognition, John Wiley & Sons, 
% New York, 2002.
%
% SEE ALSO
% MAPPINGS, DATASETS, REGOPTC, NMC, NMSC, LDC, UDC, QUADRC, NORMAL_MAP

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: qdc.m,v 1.7 2008/03/20 09:25:10 duin Exp $

function w = generalized_qdc(a, bRegularisation, alfa)

	prtrace(mfilename);

	if (nargin < 3)
		prwarning(4,'Regularisation parameter alfa not given, estimating alfa from class frequencies.');
		alfa = [];
	end
	
	if (nargin < 2) 
		prwarning(4,'Regularisation not specified, assuming NO regularisation.');
		alfa = 0;
        bRegularisation = false;
	end
	
	if (nargin < 1)       % No input arguments: 
		w = mapping(mfilename,{bRegularisation,alfa}); % return an untrained mapping.
		
	else % training
		
		islabtype(a,'crisp','soft'); % Assert A has the right labtype.
		isvaldfile(a,2,2); % at least 2 objects per class, 2 classes

		[m,k,c] = getsize(a);

%         weights = getident(a,'weights');
%         
%         if( isempty(weights) )
%             weights = ones(m,1);
%         end
        
		[U,G] = meancov(a);

		% Calculate means and priors.

		pars.mean  = +U;
		pars.prior = getprior(a);

        if(bRegularisation)         
            Nk = classsizes(a);
            [tmpU,Gldc] = meancov(seldat(a,1)); 
            Gldc = pars.prior(1)*Gldc;
            for i=2:c
                [tmpU,tmpG] = meancov(seldat(a,i)); 
                Gldc = Gldc + pars.prior(i)*tmpG;
            end
            clear tmpU tmpG;
        end
        
		% Calculate class covariance matrices.

		pars.cov   = zeros(k,k,c);
		for j = 1:c
			F = G(:,:,j);
		
			% Regularize, if requested.
			if (bRegularisation)
                
                if( isempty(alfa) || isnan(alfa) || isinf(alfa) )
                    % Usamos una exponencial para mapear el ratio m/k a un
                    % valor de alfa. Sabemos que para n/k 10 estamos en una
                    % situación adecuada para el calculo de la mat. de cov.
                    % En ese punto seteamos un porcentaje de mezcla, por
                    % ejemplo 70 de G_k y 30 de Gldc
%                     alfa_used = exp(-Nk(j)/k/1000);
                    
                    if Nk(j) < 10000
                        alfa_used = 1;
                    else
                        alfa_used = 0;
                    end
                else
                    alfa_used = alfa;
                end
                
				F = ((1-alfa_used) * Nk(j)* F + alfa_used * m * Gldc)/((1-alfa_used) * Nk(j) + alfa_used * m );
                
			end

			pars.cov(:,:,j) = F;
		end

		w = normal_map(pars,getlab(U),k,c);
		w = setcost(w,a);
		
	end

	w = setname(w,'Bayes-Normal-3');

return;
