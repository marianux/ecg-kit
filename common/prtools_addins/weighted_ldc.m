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
% covariance matrix is the weighted (by a priori probabilities) average of
% the class covariance matrices. R and S (0 <= R,S <= 1) are regularization
% parameters used for finding the covariance matrix G by:
%
%      G = (1-R-S)*G + R*diag(diag(G)) + S*mean(diag(G))*eye(size(G,1))
%
% This covariance matrix is then decomposed as G = W*W' + sigma^2 * eye(K),
% where W is a K x M matrix containing the M leading principal components
% and sigma^2 is the mean of the K-M smallest eigenvalues.
%
% The use of soft labels is supported. The classification A*W is computed
% by NORMAL_MAP.
%
% If R, S or M is NaN the regularisation parameter is optimised by REGOPTC.
% The best result are usually obtained by R = 0, S = NaN, M = [], or by
% R = 0, S = 0, M = NaN (which is for problems of moderate or low dimensionality
% faster). If no regularisation is supplied a pseudo-inverse of the
% covariance matrix is used in case it is close to singular.
%
% Note that A*(KLMS([],N)*NMC) performs a similar operation by first
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

function [W,r,s,dim] = weighted_ldc(a, r, s, dim, bUseWeight_Vector)

	prtrace(mfilename);

    if (nargin < 5) || isempty(bUseWeight_Vector)
        prwarning(4,'Weighting parameter W not given, assuming equal weights');
        bUseWeight_Vector = false;
    end
    
	if (nargin < 4) || isempty(dim)
        prwarning(4,'subspace dimensionality M not given, assuming K');
        dim = [];
    end
	if (nargin < 3) || isempty(s)
		prwarning(4,'Regularisation parameter S not given, assuming 0.');
		s = 0; 
    end
    
	if (nargin < 2) || isempty(r)
		prwarning(4,'Regularisation parameter R not given, assuming 0.');
		r = 0;
    end

	if (nargin < 1) | (isempty(a))  % No input arguments: 
		W = prmapping(mfilename,{r,s,dim,bUseWeight_Vector}); % return an untrained mapping.
	
	elseif any(isnan([r,s,dim]))    % optimize regularisation parameters
			defs = {0,0,[]};
			parmin_max = [1e-8,9.9999e-1;1e-8,9.9999e-1;1,size(a,2)];
			[W,r,s,dim] = regoptc(a,mfilename,{r,s,dim,bUseWeight_Vector},defs,[3 2 1],parmin_max,testc([],'soft'),[1 1 0]);
	
	else % training

		islabtype(a,'crisp','soft');
		isvaldfile(a,2,2); % at least 2 object per class, 2 classes

		[m,k,c] = getsize(a);
        X = (+a)';
%         not_nan_idx = find(all(~isnan(X)));
%         m  = length(not_nan_idx);

        if (nargin < 2) | isempty(bUseWeight_Vector)
            prwarning(4,'Weighting parameter W not given, assuming equal weights');
            weight_matrix = ones(c,m);
            bUseWeight_Vector = false;
        else
            
            if( bUseWeight_Vector )
                feature_vector_weight = getident(a, 'feature_vector_weight');
            end
            
            weight_matrix = zeros(c,m);
            for ii = 1:c
                class_indexes = findnlab(a,ii);
                if( bUseWeight_Vector )
                    weight_matrix(ii,class_indexes) = feature_vector_weight(class_indexes);
                else
                    weight_matrix(ii,class_indexes) = 1;
                end
            end

            if(any(isnan(weight_matrix)))
                error('Something VERY wrong with W');
            end
            
        end
        
        % Averiguo los features direccionales.
        cFeaturesDomain = getfeatdom(a);

        DirectionalFeatures = [];

        for ii = 1:length(cFeaturesDomain)
            if(~isempty(cFeaturesDomain{ii}))
                DirectionalFeatures = [DirectionalFeatures ii];
            end
        end
        
        Allfeatures_idx = 1:k;
        
        if( ~isempty(DirectionalFeatures) )
            OtherFeatures = setdiff(Allfeatures_idx,DirectionalFeatures);
        else
            OtherFeatures = Allfeatures_idx;
        end
        
        U = nan(k,c);
        
%         if( ~isempty(OtherFeatures) )
%             U(OtherFeatures,:)=X(OtherFeatures,not_nan_idx)*weight_matrix(:,not_nan_idx)'/K;
            for ii = rowvec(OtherFeatures)
                not_nan_idx = find(~isnan(X(ii,:)));
                K=diag(sum(weight_matrix(:,not_nan_idx),2));
                U(ii,:)=X(ii,not_nan_idx)*weight_matrix(:,not_nan_idx)'/K;
            end
%         end
        
%         if( ~isempty(DirectionalFeatures) )
%             U(DirectionalFeatures,:)=angle(exp(1i*X(DirectionalFeatures,not_nan_idx))*weight_matrix(:,not_nan_idx)'/K);
            for ii = rowvec(DirectionalFeatures)
                not_nan_idx = find(~isnan(X(ii,:)));
                K=diag(sum(weight_matrix(:,not_nan_idx),2));
                U(ii,:)=angle(exp(1i*X(ii,not_nan_idx))*weight_matrix(:,not_nan_idx)'/K);
            end
            
%         end
        
        
%         U = nan(k,c);
%         for ii = 1:c
%             class_indexes = findnlab(a,ii);
%             U(:,ii)=median(X(:,class_indexes),2);
%         end

        iAux = nan(k, m);
        not_nan_idx = 1:m;
        weight_matrix_sqrt = sqrt(weight_matrix(:,not_nan_idx));
        K=diag(sum(weight_matrix,2));

        if( ~isempty(OtherFeatures) )
            iAux(OtherFeatures,:) = (  ( X(OtherFeatures,not_nan_idx).*(ones(length(OtherFeatures),c)*weight_matrix_sqrt) ) - U(OtherFeatures,:)*weight_matrix_sqrt );
        end

        if( ~isempty(DirectionalFeatures) )        
            for ii = DirectionalFeatures                                     
                iAux1 = ( X(ii,:).*(ones(1,c)*weight_matrix_sqrt) );
                iAux2 = U(ii,:)*weight_matrix_sqrt;
                iAux3 = [ iAux1 - iAux2; 2*pi + iAux1 - iAux2; iAux1 - 2*pi - iAux2 ];
                [dummy, iMinIndex] = min(abs(iAux3), [], 1);
                iMinIndex = sub2ind(size(iAux3),iMinIndex,1:m);
                iAux(ii,:) = iAux3(iMinIndex);
            end 
        end
        
%         G = iAux * iAux' / sum(sum(K));
        G = zeros(k,k);
        for ii = 1:k
            for jj = 1:k
                aux_val = iAux(ii,:) .* iAux(jj,:);
                bAux = ~isnan(aux_val);
                if( ~any(bAux) ); continue; end;
                K = sum(sum(weight_matrix(:,bAux)));
                G(ii,jj) = sum(aux_val(bAux)) / K;
            end
        end
        
        %Cuando no separaba features direccionales de no direccionales.
%         G=X*(  ( X.*(ones(k,c)*weight_matrix) ) - U*weight_matrix  )' / sum(sum(K));
    
		% Regularize 
		if (s > 0) || (r > 0)
			G = (1-r-s)*G + r * diag(diag(G)) + s*mean(diag(G))*eye(size(G,1));
		end	
	
        if (dim < k)
                dim = min(rank(G)-1,dim);
            [eigvec,eigval] = eig(G); eigval = diag(eigval);
            [dummy,ind] = sort(-eigval);
                % Estimate sigma^2 as avg. eigenvalue outside subspace.
            sigma2 = mean(eigval(ind(dim+1:end)));
            % Subspace basis: first DIM eigenvectors * sqrt(eigenvalues).
            G = eigvec(:,ind(1:dim)) * diag(eigval(ind(1:dim))) * eigvec(:,ind(1:dim))' + ...
              sigma2 * eye(k);
        end

        w.mean  = U';
		w.prior = getprior(a);
		w.cov = G;
		W = normal_map_new(w,getlablist(a),k,c);
		W = setcost(W,a);
        W = setuser(W,cFeaturesDomain);
		
	end
	
	W = setname(W,'Bayes-Normal-1');
	
return
