%PARZEN_MAP Map a dataset on a Parzen densities based classifier
% 
% 	F = PARZEN_MAP(A,W)
%
% INPUT
%   A   Dataset
%   W   Trained Parzen classifier mapping (default: PARZENC(A))
%
% OUTPUT
%   F   Mapped dataset
%
% DESCRIPTION 
% Maps the dataset A by the Parzen density based classfier W. F*sigm are the
% posterior probabilities. W should be trained by a classifier like PARZENC.
% This routine is called automatically to solve A*W, if W is trained by
% PARZENC.
% 
% The global PRMEMORY is read for the maximum size of the internally declared 
% matrix, default inf.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, PARZENC, TESTP

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: parzen_map.m,v 1.5 2008/08/18 21:16:26 duin Exp $

function f = parzen_map(b,w)

		% If no mapping is supplied, train one.
	if (nargin < 2)
		w = parzenc(b); 
	end

	pars = getdata(w);
	a = pars{1}; 								% Stored training dataset.
	h = pars{2}; 								% Smoothing paramater.

	[m,k,c] = getsize(a); nlab = getnlab(a); p = getprior(a);

	% do we have priors set in the mapping?
	if length(pars)==3
	   p=pars{3};
	end
	
	% Prepare a matrix with smoothing parameters for each class and feature.

	if (size(h,1) == 1)
		h = repmat(h,c,1);
	end
	if (size(h,2) == 1)
		h = repmat(h,1,k);
	end
	
	
	if (any(size(h) ~= [c,k]))
		error('The size of the array with smoothing parameters does not match that of the training set of W.');
	end

	[mt,kt] = size(b);
	if (kt ~= k)
		error('The size of the set A does not match that of the training set of W.'); 
	end

	[num,n] = prmem(mt,m);
	f = ones(mt,c); 					% Prob. densities for each test sample and class.

	for j = 0:num-1						% Process each chunk of the test set separately.

		if (j == num-1)
			nn = mt - num*n + n;	% Last chunk may have smaller size.
		else
			nn = n;
		end
		range = [j*n+1:j*n+nn];

		% Estimate the class probabilities.

		for i = 1:c
			
			if islabtype(a,'crisp')
				I = findnlab(a,i); hh = h(i,:);		% Parameters for this class.

				% Calculate squared distances to kernel centers.
				
				D = +distm(a(I,:)./repmat(hh,length(I),1), ...
				           +b(range,:)./repmat(hh,length(range),1));
				if (length(I) > 0)
					f(range,i) = mean(exp(-D*0.5),1)';			% Apply kernel.
				end
				const = repmat(-log(length(I)),length(I),1);
			else
				hh = h(i,:);
				I = find(a.targets(:,i) > 0); % Avoid objects with zero weight
				v = a.targets(I,i);
				D = distm(+a(I,:)./repmat(hh,length(I),1), ...
                 +b(range,:)./repmat(hh,length(range),1));
				f(range,i) = exp(-D'*0.5) * v / sum(v);
				const = log(v/sum(v));
			end
			% Normalize and multiply by class prior to get a density. Add REALMIN
			% to prevent division-by-zero errors.
			f(range,i) = p(i)*f(range,i)/((sqrt(2*pi).^k)*prod(hh)+realmin);
			if (getout_conv(w) == 2) % take log of density to preserve tails
				J = find(f(range,i) <= 1e-303);
				N = find(f(range,i) >  1e-303);
				f(range(N),i) = log(f(range(N),i));
				[dm,R] = min(D(:,J)); % for zero densities use nearest neighbor only
				%f(range(J),i) = log(p(i)/((sqrt(2*pi).^k)*prod(hh)+realmin)) - dm'*0.5 + const(R);
				f(range(J),i) = log(p(i)) - k*log(sqrt(2*pi)) - sum(log(hh)) - dm'*0.5 + const(R);
				% in case of soft labels this is tricki. We just hope that the
				% weight of the nearest neighbor is sufficiently large
			end
		end
	end
	if (getout_conv(w) == 2) % scale to gain accuracy in the tails
		fmax = max(f,[],2);
		f = f - repmat(fmax,1,c);
		f = exp(f);
	else
		f = f + realmin; % avoid devision by 0 in computing posterios later
  end
	
  if isdataset(b)
    f = setdata(b,f,getlabels(w));
  end
return
