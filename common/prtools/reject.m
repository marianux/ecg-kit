%REJECT Compute the error-reject trade-off curve
% 
%   E = REJECT(D);
%   E = REJECT(A,W);
%	
% INPUT
%   D   Classification result, D = A*W
%   A   Dataset
%   W   Cell array of trained classifiers
%
% OUTPUT
%   E   Structure storing the error curve and information needed for plotting
%
% DESCRIPTION
% E = REJECT(D) computes the error-reject curve of the classification 
% result D = A*W, in which A is a dataset and W is a classifier. E is 
% a structure storing the error curve in E.ERROR. Use PLOTE(E) for 
% plotting the result.
%
%	E = REJECT(A,W) computes a set of error-reject curves for all trained 
% classifiers stored in the cell array W.
%
% EXAMPLES
% A - training set, B - test set:
% D = B*NMC(A); E = REJECT(D); PLOTE(E);   % Plots a single curve
% E = REJECT(B,A*{NMC,UDC,QDC}); PLOTE(E); % Plots 3 curves
% 
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd edition, 
% John Wiley and Sons, New York, 2001.
% 2. A. Webb, Statistical Pattern Recognition, John Wiley & Sons, New York, 2002.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PLOTE, PRROC, TESTC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: reject.m,v 1.3 2007/02/15 10:11:15 davidt Exp $

function e = reject(A,W)
		% No data, return an untrained mapping.
	if (nargin == 0) | (isempty(A))
		e = prmapping('reject','fixed');
		return;
	end

	compute = 0;					% COMPUTE is 1 if the classier W should be applied to A.
	name = [];
	if (nargin == 2)
		if iscell(W)
			classfno = length(W); 	% Number of classifiers, hence error curves.
			compute = 1;
		elseif ismapping(W) % W is one classifier only.
			if ~istrained(W)
				error('Classifier should be trained.')
			end
			D = A*W*classc;   % Normalize the outputs of W and apply it to the data A.
			classfno = 1;			
			name = getname(W);
		else
			error('Second parameter should be a classifier or a cell array of classifiers.')
		end
	else									% Only one input argument. 
		D = A;
		classfno = 1; 			
	end

	% Fill the structure E with information needed for plotting.
	m = size(A,1);
	e.error = zeros(classfno,m+1);
  % m+1 for storing error of 0 size dataset
	e.std   = zeros(classfno,m+1);
	e.xvalues = [0:m]/m;;
	e.n = 1;
	datname = getname(A);
	if ~isempty(datname)
		e.title = ['Reject curve for the ' datname];
	end
	e.xlabel= 'Reject';
	e.ylabel= 'Error';
	e.plot  = 'plot';
	e.names = name;

	for j=1:classfno
		if (compute)
			w = W{j};
			if ~istrained(w)
				error('Classifier should be trained.')
			end
			D = A*w*classc; 		% Normalize the outputs of W and apply it to A.
			e.names = char(e.names,getname(w));
		end
		% Compare the classification result with the actual labels.
		% N is a 0/1 vector pointing to all distinct/equal labels.
		[err,n] = nlabcmp(labeld(D),getlab(D)); 

		% A trick: if D consists of one column, as before returned by 
		% FISHERC in case of a 2-class problem. 
		% May be they don't exist anymore in PRTools, but once they did
		% and may be they pop up again.
		if (size(D,2) == 1)
			D = [D 1-D];
		end
		% Sort the objects, starting from the closest to the decision 
		% boundary to the furthest away.
		[y,J] = sort(max(+D,[],2));

		% 1/0 in N corresponds now to objects correctly/wrongly classified.
		n = 1-n(J)';				 
		e.error(j,:) = [err err-cumsum(n)]/m;
	end
	if (classfno > 1)
		e.names(1,:) = [];
	end

return;
