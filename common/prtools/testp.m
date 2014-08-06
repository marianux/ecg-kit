%TESTP Error estimation of Parzen classifier
% 
% 	E = TESTP(A,H,T)
% 	E = TESTP(A,H)
% 
% INPUT
%     A    input dataset
%     H    matrix smoothing parameters (optional, def: determined via
%            parzenc)
%     T    test dataset (optional)
%
% OUTPUT
%     E    estimated error rate
%
% DESCRIPTION 
% Tests a dataset T on dataset A using a Parzen classification and returns
% the classification error E.  Returns the leave-one-out error estimate. If
% H is not given, it is determined by PARZENC.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PARZEN_MAP, PARZENML, PARZENC. 

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% may be to be merged with parzen_map, see also testk

% $Id: testp.m,v 1.4 2009/09/08 21:27:51 duin Exp $

function [e,d] = testp(a,h,t)

		
	isvaldfile(a,2,2);
	a = testdatasize(a);
	a = testdatasize(a,'objects');
	
	if nargin < 2, [W,h] = parzenc(a); end
	[m,k,c] = getsize(a);
	nlab = getnlab(a);
	lablist = getlablist(a);
	p = getprior(a);
	
	if length(h) == 1, h = h*ones(1,c); end
	if length(h) ~= c, error('Wrong number of smoothing parameters'); end

	% if no test dataset is specified
	if nargin <= 2
		% find for each sample cross-validated estimate
		% of aposteriori probability.
		d = classp(a,nlab,h,p);
		[dmax,J] = max(d',[],1);
		e = nlabcmp(lablist(J,:),lablist(nlab,:)) / m;
        % if the validation dataset is given
	elseif nargin == 3
		lablistt = getlablist(t);
		[n,kt] = size(t);
		nlabt = getnlab(t);
		if k ~= kt 
			error('Data sizes do not match');
		end
		d = classp(a,nlab,h,p,t); 
		[dmax,J] = max(d',[],1);
		e = nlabcmp(lablistt(J,:),lablistt(nlabt,:)) / n;
	end
return

%CLASSP estimate of Parzen density (if t is not specified, 
% a leave-one-out error estimate for a is returned)
function F = classp(a,nlab,h,p,t)

	[m,k] = size(a);
	maxa = max(max(abs(a)));
	a = a/maxa;
	h = h/maxa;
	if nargin < 5
		mt = m;
	else
		[mt,kt] = size(t);
		t = t/maxa;
	end

	c = max(nlab);
	alf=sqrt(2*pi)^k; % density normalization factor
	[num,n] = prmem(mt,m); % use batches to avoid excessive memory usage
	F = ones(mt,c);
  s = sprintf('Compute distance matrix in %i batches: ',num);
  prwaitbar(num,s);
	for i = 0:num-1
    prwaitbar(num,i+1,[s int2str(i+1)]);				
		if i == num-1
			nn = mt - num*n + n;
		else
			nn = n;
		end
		range = [i*n+1:i*n+nn];
		if nargin <= 4
			D = +distm(a,a(range,:));
			D(i*n+1:m+1:i*n+nn*m) = inf*ones(1,nn); % set distances to itself at inf
		else
			D = +distm(a,t(range,:));
		end
		for i=1:c
			I = find(nlab == i);
			if length(I) > 0
				F(range,i) = p(i).*sum(exp(-D(I,:)*0.5./(h(i).^2)),1)'./(length(I)*alf*h(i)^k);
			end
		end
  end
  prwaitbar(0);
	F = F + realmin;
	F = F ./ (sum(F')'*ones(1,c));

return
