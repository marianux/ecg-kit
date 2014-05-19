%MEANCOV Estimation of the means and covariances from multiclass data
% 
%   [U,G] = meancov_new(A,N)
% 
%	INPUT
% 	A		Dataset
%   N		Normalization to use for calculating covariances: by M, the number
%					of samples in A (N = 1) or by M-1 (default, unbiased, N = 0).
%
% OUTPUT
% 	U		Mean vectors
% 	G		Covariance matrices
%
% DESCRIPTION  
%	Computation of a set of mean vectors U and a set of covariance matrices G
%	of the C classes in the dataset A. The covariance matrices are stored as a
%	3-dimensional matrix G of the size K x K x C, the class mean vectors as a
%	labeled dataset U of the size C x K.
%
% The use of soft labels or target labels is supported.
% 
%	SEE ALSO 
%	DATASETS, GAUSS, NBAYESC, DISTMAHA

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: meancov.m,v 1.17 2005/10/19 08:32:11 duin Exp $

function [U,G] = meancov_new(a,n)

	prtrace(mfilename);

	% N determines whether the covariances are normalized by M (N = 1) or by 
	% M-1 (unbiased, N = 0), where M is the number of objects.

	if (nargin < 2)
		prwarning(4,'normalisation not specified, assuming by M-1');
		n = 0;
	end

	if (n ~= 1) && (n ~= 0)
		error('Second parameter should be either 0 or 1.')
    end

	if (~isa(a,'prdataset'))			% A is a matrix: compute mean and covariances
        
        DirectionalFeatures = [];
        
        if( isempty(DirectionalFeatures) )
            U = nanmean(a);
            G = nancov(a,n);
        else
            
            iCantFeatures = size(a,2);
            U = nan(1, iCantFeatures);
            NotDirectionalFeatures = setdiff(1:iCantFeatures, DirectionalFeatures);

            if( ~isempty(NotDirectionalFeatures) )
                U(NotDirectionalFeatures) = nanmean(a(:,NotDirectionalFeatures));							% 	in the usual way.
            end
            
            U(DirectionalFeatures) = angle(nanmean(exp(1i*a(:,DirectionalFeatures))));
            
            bNotNaN = ~any(isnan( a(:, DirectionalFeatures ) ),2);
            iCantElementos = sum(bNotNaN);
            
            if(iCantElementos < 2)
                
                G = zeros(iCantFeatures);
                
            else
                iAux = nan(iCantElementos, iCantFeatures);

                if( ~isempty(NotDirectionalFeatures) )
                    iAux(:,NotDirectionalFeatures) = a(bNotNaN,NotDirectionalFeatures) - repmat(U(NotDirectionalFeatures), iCantElementos,1);
                end

                for jj = DirectionalFeatures                                     
                    iAux1 = a(bNotNaN, jj );
                    iAux2 = repmat(U(jj), iCantElementos,1);
                    iAux3 = [ iAux1 - iAux2, 2*pi + iAux1 - iAux2, iAux1 - 2*pi - iAux2 ];
                    [dummy, iMinIndex] = min(abs(iAux3), [], 2);
                    iMinIndex = sub2ind(size(iAux3),1:size(iAux3,1),iMinIndex');
                    iAux(:,jj) = iAux3(iMinIndex);
                end 

                if(n == 1)
                    G = 1/iCantElementos*(iAux'*iAux);
                else
                    G = 1/(iCantElementos-1)*(iAux'*iAux);
                end
            end

        end
        
    else
        
		[m,k,c] = getsize(a);

        cFeaturesDomain = getfeatdom(a);
        
        DirectionalFeatures = [];
        
        for jj = 1:length(cFeaturesDomain)
            if(~isempty(cFeaturesDomain{jj}))
                DirectionalFeatures = [DirectionalFeatures jj];
            end
        end

        if (islabtype(a,'crisp'))
			
			if (c==0)
                
                aa = +a;
                
                if( isempty(DirectionalFeatures) )
                    U = nanmean(aa);
                    G = nancov(aa,n);
                else

                    iCantFeatures = k;
                    U = nan(1, iCantFeatures);
                    NotDirectionalFeatures = setdiff(1:iCantFeatures, DirectionalFeatures);

                    if( ~isempty(NotDirectionalFeatures) )
                        U(NotDirectionalFeatures) = nanmean(aa(:,NotDirectionalFeatures));							% 	in the usual way.
                    end
                    
                    U(DirectionalFeatures) = angle(nanmean(exp(1i*aa(:,DirectionalFeatures))));

                    bNotNaN = ~any(isnan( aa(:, DirectionalFeatures ) ),2);
                    iCantElementos = sum(bNotNaN);
                    iAux = nan(iCantElementos, iCantFeatures);
                    
                    if( ~isempty(NotDirectionalFeatures) )
                        iAux(:,NotDirectionalFeatures) = aa(bNotNaN,NotDirectionalFeatures) - repmat(U(NotDirectionalFeatures), iCantElementos,1);
                    end


                    for jj = DirectionalFeatures                                     
                        iAux1 = aa(bNotNaN, jj );
                        iAux2 = repmat(U(jj), iCantElementos,1);
                        iAux3 = [ iAux1 - iAux2, 2*pi + iAux1 - iAux2, iAux1 - 2*pi - iAux2 ];
                        [dummy, iMinIndex] = min(abs(iAux3), [], 2);
                        iMinIndex = sub2ind(size(iAux3),1:size(iAux3,1),iMinIndex');
                        iAux(:,jj) = iAux3(iMinIndex);
                    end 

%                     iAux1 = aa(:,DirectionalFeatures) - repmat(U(DirectionalFeatures), iCantElementos,1);
%                     boolAbs = abs( iAux1 ) > pi;
%                     SignX = sign( iAux1 );
% 
%                     bAux = boolAbs & SignX < 0;
% 
%                     aa( bAux, DirectionalFeatures ) = aa( bAux, DirectionalFeatures ) + 2*pi;
%                     aa( bAux, DirectionalFeatures ) = aa( bAux, DirectionalFeatures ) - 2*pi;
% 
%                     iAux(:,DirectionalFeatures) = rem( aa(:, DirectionalFeatures ) - repmat(U(DirectionalFeatures), iCantElementos,1) , 2*pi);

                    if(n == 1)
                        G = 1/iCantElementos*(iAux'*iAux);
                    else
                        G = 1/(iCantElementos-1)*(iAux'*iAux);
                    end

                end
                
            else

                U = nan(c,k);
                G = nan(k,k,c);
                
                if( isempty(DirectionalFeatures) )
                    
                    for ii = 1:c     
                        J = findnlab(a,ii);
                        if( ~isempty(J) )
                            aa = a(J,:);
                            bNotNaN = ~any(isnan( aa ),2);
                            U(ii,:) = nanmean( aa(bNotNaN,:) ,1);
                            if (nargout > 1 && sum(bNotNaN) > 1 )
                                G(:,:,ii) = covm(aa(bNotNaN,:),n);	
                            else
                                G(:,:,ii) = 0;
                            end
                        end
                    end
                    
                else
                
                    NotDirectionalFeatures = setdiff(1:k, DirectionalFeatures);
                    aa = +a;

                    U = nan(c, k);
                    
                    for ii = 1:c     
                        J = findnlab(a,ii);
                        if( ~isempty(J) )

                            if( ~isempty(NotDirectionalFeatures) )
                                U(ii,NotDirectionalFeatures) = nanmean(aa(J,NotDirectionalFeatures));							% 	in the usual way.
                            end
                            
                            U(ii,DirectionalFeatures) = angle(nanmean(exp(1i*aa(J,DirectionalFeatures))));
                            
                            bNotNaN = ~any(isnan( aa(J,:) ),2);
                            iCantElementos = sum(bNotNaN);
                            
                            if(iCantElementos < 2)
                                G(:,:,ii) = zeros(k);
                            else
                                
                                if (nargout > 1)

                                    iAux = nan(iCantElementos, k);

                                    if( ~isempty(NotDirectionalFeatures) )
                                        iAux(:,NotDirectionalFeatures) = aa(J(bNotNaN),NotDirectionalFeatures) - repmat(U(ii,NotDirectionalFeatures,:), iCantElementos,1);
                                    end

                                    for jj = DirectionalFeatures
                                        iAux1 = aa(J(bNotNaN), jj );
                                        iAux2 = repmat(U(ii,jj), iCantElementos,1);
                                        iAux3 = [ iAux1 - iAux2, 2*pi + iAux1 - iAux2, iAux1 - 2*pi - iAux2 ];
                                        [dummy, iMinIndex] = min(abs(iAux3), [], 2);
                                        iMinIndex = sub2ind(size(iAux3),1:size(iAux3,1),iMinIndex');
                                        iAux(:,jj) = iAux3(iMinIndex);
                                    end

                                    if(n == 1)
                                        G(:,:,ii) = 1/iCantElementos*(iAux'*iAux);
                                    else
                                        G(:,:,ii) = 1/(iCantElementos-1)*(iAux'*iAux);
                                    end
                                end
                            end
                        end
                    end
                end
			end
			labu = getlablist(a);
  	elseif (islabtype(a,'soft'))
  		problab = gettargets(a);
			% Here we also have to be careful for unlabeled data
			if (c==0)
				prwarning(2,'The dataset has soft labels but no targets defined: using targets 1');
				U = nanmean(+a);
				G = nancov(+a,n);
			else
				U = zeros(c,k);
				for ii = 1:c

					% Calculate relative weights for the means.
					g = problab(:,ii); nn = sum(g); g = g/mean(g); 

					U(ii,:) = mean(a.*repmat(g,1,k));	% Weighted mean vectors	

					if (nargout > 1)

						u  = mean(a.*repmat(sqrt(g),1,k));

						% this appears to be needed to weight cov terms properly
						G(:,:,ii) = covm(a.*repmat(sqrt(g),1,k),1) - U(ii,:)'*U(ii,:) + u'*u;

						% Re-normalise by M-1 if requested.
						if (n == 0)
							G(:,:,ii) = m*G(:,:,ii)/(m-1);
						end
					end
				end
            end
			labu = getlablist(a);
        else
			% Default action.
            U = mean(a);
            G = covm(a,n);
			labu = [];
        end

		% Add attributes of A to U.
		U = prdataset(U,labu,'featlab',getfeatlab(a), ...
                                'featsize',getfeatsize(a));
		if (~islabtype(a,'targets'))
			p = getprior(a);
			U = setprior(U,p); 
		end

	end;

return
