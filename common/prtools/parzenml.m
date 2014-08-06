%PARZENML Optimum smoothing parameter in Parzen density estimation.
% 
%   H = PARZENML(A)
% 
% INPUT	
%   A    Input dataset
%
% OUTPUT
%   H    Scalar smoothing parameter (in case of crisp labels)
%        Vector with smoothing parameters (in case of soft labels)
%
% DESCRIPTION
% Maximum likelihood estimation for the smoothing parameter H in the 
% Parzen denstity estimation of the data in A. A leave-one out 
% maximum likelihood estimation is used. 
%
% The dataset A can either be crisp or soft labeled. In case of crisp
% labeling the class information is not used and a single smoothing 
% parameter is estimated. In case of soft labels a smoothing parameter
% for every class is estimated and objects are weighted in relation to
% their class weigthts (soft label value). 
% It may be profitable to scale the data before calling it. eg. 
% WS = SCALEM(A,'variance'); A = A*WS.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, SCALEM, SELDAT, PARZENM, PARZENDC, PRPROGRESS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: parzenml.m,v 1.11 2010/03/25 15:39:46 duin Exp $

function h = parzenml(A,fid)

		
	if nargin < 2, fid = []; end

	if isdouble(A), A = prdataset(A); end
	
	A = testdatasize(A);
	A = testdatasize(A,'objects');

	if islabtype(A,'crisp')
		h = parzenmlc(A,fid);
	elseif islabtype(A,'soft')
		h = parzenmls(A,fid);
	else
		error('Label type should be either ''crisp'' or ''soft''')
	end
	
	return
	
function h = parzenmlc(A,fid) %crisp version

	[m,k] = size(A);
	DD= distm(+A) + diag(1e70*ones(1,m));
	E = min(DD);
	
	h1 = sqrt(max(E));    % initial estimate of h
	F1 = derlc(DD,E,h1,k); % derivative

	prprogress(fid,'parzenml:\n');
	prprogress(fid,' %6.4f   %6.3e\n',h1,F1);
	if abs(F1) < 1e-70 
		h = h1;
		prwarning(4,'jump out\n');
		return;
	end
	
	a1 = (F1+m*k)*h1*h1;
	h2 = sqrt(a1/(m*k));  % second guess
	F2 = derlc(DD,E,h2,k); % derivative

	prprogress(fid,' %6.4f   %6.3e\n',h2,F2);
	if (abs(F2) < 1e-70) | (abs(1e0-h1/h2) < 1e-6) 
		h = h2;
		prwarning(4,'jump out\n');
		return
	end
	
	% find zero-point of derivative to optimize h^2
	% stop if improvement is small, or h does not change significantly
	
	alf = 1;
	prwaitbar(100,'parzenml: Optimizing smoothing parameter',m > 100)
	iter = 0;
	while abs(1e0-F2/F1) > 1e-4 & abs(1e0-h2/h1) > 1e-3 & abs(F2) > 1e-70
		iter = iter+1;
		h3 = (h1*h1*h2*h2)*(F2-F1)/(F2*h2*h2-F1*h1*h1);
		if h3 < 0 % this should not happen
			h3 = sqrt((F2+m*k)*h2*h2/(m*k));
		else
			h3 = sqrt(h3);
		end
		prwaitbar(100,100-100*exp(-iter/10));
		h3 = h2 +alf*(h3-h2);
		F3 = derlc(DD,E,h3,k);
		prprogress(fid,' %6.4f   %6.3e\n',h3,F3);
		F1 = F2; F2 = F3;
		h1 = h2; h2 = h3;
		alf = alf*0.99; % decrease step size
	end
	h = h2;
	prwaitbar(0);

return

function F = derlc(DD,E,h,k) % crisp version
	% computation of the likelihood derivative for Parzen density
	% given distances D and their object minima E (for increased accuracy)
	m = size(DD,1);
	warning off MATLAB:divideByZero;
		Y = (DD-repmat(E,m,1))/(2*h*h); % correct for minimum distance to save accuracy
	warning on MATLAB:divideByZero;
	IY = find(Y<20);                % take small distance only, others don't contribute
	P = zeros(m,m);
	P(IY) = exp(-Y(IY));
	PP = sum(P,2)';
	FU = repmat(realmax,1,m);
	J = find(PP~=0); 
	FU(J) = 1./PP(J);
	FF = sum(DD.*P,2);
	warning off MATLAB:divideByZero;
		F = (FU*FF)./(h*h) - m*k;
	warning on MATLAB:divideByZero;
return

function h = parzenmls(A,fid) %soft version

	SS = gettargets(setlabtype(A,'soft'));
	[m,k,c] = getsize(A);
	DD= distm(+A) + diag(1e70*ones(1,m));
	E = min(DD);
	h = zeros(c,1);
	h0 = sqrt(max(E));    % initial estimate of h
	
	
	s = sprintf('parzenml: runover classes');
	prwaitbar(c,s,m > 100);
	iter = 0;
	
	for j=1:c
		prwaitbar(c,j)
		S = SS(:,j);
		h1 = h0;
		F1 = derls(DD,E,h1,k,S); % derivative

	  prprogress(fid,'parzenml: class %i : \n',j);
		prprogress(fid,' %6.4f   %6.3e\n',h1,F1);
		if abs(F1) < 1e-70 
			h(j) = h1;
			prwarning(4,'jump out\n');
			break;
		end
	
		a1 = (F1+m*k)*h1*h1;
		h2 = sqrt(a1/(m*k));  % second guess
		F2 = derls(DD,E,h2,k,S); % derivative

		prprogress(fid,' %6.4f   %6.3e\n',h2,F2);
		if (abs(F2) < 1e-70) | (abs(1e0-h1/h2) < 1e-6) 
			h(j) = h2;
			prwarning(4,'jump out\n');
			break;
		end
	
		% find zero-point of derivative to optimize h^2
		% stop if improvement is small, or h does not change significantly
	
		
		prwaitbar(100,'parzenml: Optimizing smoothing parameter',m > 100)
		iter = 0;
		alf = 1;
		while abs(1e0-F2/F1) > 1e-4 & abs(1e0-h2/h1) > 1e-3 & abs(F2) > 1e-70
			iter = iter+1;
			prwaitbar(100,100-100*exp(-iter/10));
			h3 = (h1*h1*h2*h2)*(F2-F1)/(F2*h2*h2-F1*h1*h1);
			if h3 < 0 % this should not happen
				h3 = sqrt((F2+m*k)*h2*h2/(m*k));
			else
				h3 = sqrt(h3);
			end
			h3 = h2 +alf*(h3-h2);
			F3 = derls(DD,E,h3,k,S);
			prprogress(fid,' %6.4f   %6.3e\n',h3,F3);
			F1 = F2; F2 = F3;
			h1 = h2; h2 = h3;
			alf = alf*0.99; % decrease step size
		end
		prwaitbar(0)
		h(j) = h2;
	end
  prwaitbar(0)
return

function F = derls(DD,E,h,k,S) %soft version
	% computation of the likelihood derivative for Parzen density
	% given distances D and their object minima E (for increased accuracy)
	% S are the object weigths
	c = size(S,2);                  % number of classes
	m = size(DD,1);
	Y = (DD-repmat(E,m,1))/(2*h*h); % correct for minimum distance to save accuracy
	IY = find(Y<20);                % take small distance only, others don't contribute
	F = 0;
	for j=1:c
		P = zeros(m,m);
		P(IY) = exp(-Y(IY));
		PP = S(:,j)'*P';
		FU = repmat(realmax,1,m);
		J = find(PP~=0);  
		FU(J) = S(J,j)'./PP(J);
		K = find(S(:,j)==0);
		FU(K) = zeros(1,length(K));
		FF = (DD.*P)*S(:,j);
		F = F + (FU*FF)./(h*h);
	end
	F = F - sum(S(:))*k;
return
