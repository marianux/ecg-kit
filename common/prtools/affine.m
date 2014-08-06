%AFFINE Construct affine (linear) mapping from parameters
%
%   W = AFFINE(R,OFFSET,LABLIST_IN,LABLIST_OUT,SIZE_IN,SIZE_OUT)
%   W = AFFINE(R,OFFSET,A)
%   W = AFFINE(W1,W2) 
%
% INPUT
%   R            Matrix of a linear mapping from a K- to an L-dimensional space
%   OFFSET       Shift applied after R; a row vector of the length L
%                (optional; default: zeros(1,L))
%   LABLIST_IN   Labels of the features of the input space
%                (optional; default: (1:K)')
%   LABLIST_OUT  Labels of the features of the output space, e.g. class names 
%                for linear classifiers (optional; default: (1:L)')
%   SIZE_IN      If based on images: size vector of the input dimensionality 
%                (optional; default: K)
%   SIZE_OUT     If based on images: size vector of the output  dimensionality 
%                (optional; default: L)
%   A            Dataset (LAB_IN_LIST and SIZE_IN are derived from A)
%   W1,W2        Affine mappings
%
% OUTPUT
%   W            Affine mapping
%
% DESCRIPTION  
% This is a low level basic PRTools routine, not intended for direct use.
% It efines a mapping W based on a linear transformation R and an offset. 
% R should be a [K x L] matrix describing a linear transformation from 
% a K-dimensional space to an L-dimensional space. If K=1, then R is 
% interpreted as the diagonal of an [L x L] diagonal matrix. OFFSET is 
% a row vector of the length L, added afterwards.
%
% Affine mappings are treated by PRTools in a special way. A scaling
% defined for an affine mapping, e.g. by W = SETSCALE(W,SCALE) is directly
% executed by a multiplication of the coefficients. Also, the product of 
% two affine mappings is directly converted to a new affine mapping. 
% This routine also executes W = AFFINE(W1,W2), if W1 and W2 are affine. 
% B = AFFINE(A,W), if A is a dataset and W is an affine mapping.
% Finally, the transpose of an affine mapping exists and is defined as 
% an another affine mapping.
%
% An [M x K] dataset A can be mapped as D = A*W. The result is equivalent
% to [+A, ones(M,1)]*[R; OFFSET]. The dataset D has feature labels stored 
% in LABLIST. The number of this labels should, thereby, be at least L.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: affine.m,v 1.10 2009/03/17 10:03:51 duin Exp $

function w = affine(R,offset,lablist_in, lablist_out,size_in,size_out)

  if (nargin == 1) | (~isa(offset,'prmapping'))

		% Definition of an affine mapping
		[m,k] = size(R);
		if (nargin < 6)
			prwarning(5,'SIZE_OUT is not specified. The number of columns of R, %d, is assumed.', k); 
			size_out = k; 
		end
		if (nargin < 5)
			prwarning(5,'SIZE_IN is not specified. The number of rows of R, %d, is assumed.', m); 
			size_in = m; 
		end
		if (nargin < 4)
			prwarning(5,'LABLIST_OUT is not specified, [1:%d]'' assumed.', k); 
			lablist_out = []; 
		end
		if (nargin < 3) 
			prwarning(5,'LABLIST_IN is not specified, [1:%d]'' assumed.', m); 
			lablist_in = []; 
		end
		if (nargin < 2) | (isempty(offset)) 
			prwarning(3,'OFFSET not specified, a zero vector assumed.'); 
			offset = zeros(1,k); 
		end

		% Check consistencies
		if (~isa(R,'double'))
			error('No proper transformation matrix stored.')
		end
		if (size_in == 1) & nargin < 3 										% R is a scaling vector
			size_in = size_out; 
		end	
		
		if (isempty(lablist_in))
			lablist_in = genlab(1,[1:size_in]');
		end
		
		cost = [];
		if (isa(lablist_in,'prdataset'))						% Copy labels from dataset/datafile
			cost = lablist_in.cost;
			size_in = getfeatsize(lablist_in);
			lablist_in = getfeatlab(lablist_in);
			%if isempty(lablist_in)
			%	lablist_in = num2str([1:size_in]');
			%end
			% size_out = k; % Wrong for classifiers defined for 1D datasets
		end
		
		if ~isempty(lablist_in) & (size(lablist_in,1) < m)
			error('Wrong number of input labels supplied.')
		end

		if isempty(lablist_out)
			lablist_out = genlab(1,[1:size_out]');
		end
		
		if (size(lablist_out,1) < k)
			error('Wrong number of output labels supplied.')
		end
		if any(size(offset) ~= [1,k])
			error('Offset is not a row vector of the correct size.')
		end
		
		% Store the results:
		d.rot = R;
		d.offset = offset;
		d.lablist_in = lablist_in;
		w = prmapping(mfilename,'trained',d,lablist_out,size_in,size_out);
		w = setcost(w,cost);
	
	elseif isa(R,'prmapping') 

		% Two mappings, stored in R and OFFSET, should be combined.
		w1 = R;
		w2 = offset;
		if (~isclassifier(w1)) & (~isclassifier(w2)) & (strcmp(getmapping_file(w1),'affine')) & (strcmp(getmapping_file(w2),'affine'))
			% Combine two affine mappings
			% If d1.rot or d2.rot are vectors, they have to be interpreted as
			% the diagonal matrices, unless the inner dimension does not fit.
			d1 = +w1; 
			d2 = +w2;
			if (size(d1.rot,1) == 1)				% d1.rot is a vector
				if (size(d2.rot) == 1)				% d2.rot is a vector
					d.rot = d1.rot.*d2.rot;
					d.offset = d1.offset.*d2.rot + d2.offset;
				else			  									% d2.rot is a matrix
					d.rot = repmat(d1.rot',1,size(d2.rot,2)).*d2.rot;
					d.offset = d1.offset*d2.rot + d2.offset;
				end
			else														% d1.rot is a matrix
%RD Here comes a bug fix that I needed to continue, I am not sure it
%RD is sufficient It may even introduce new problems, especially for
%   1D datasets.
				%if size(d2.rot,1) == 1 		% d2.rot is vector
				if (size(d1.rot,2) > 1) & (size(d2.rot,1) == 1) 	% d2.rot is a vector
					d.rot = d1.rot.*repmat(d2.rot,size(d1.rot,1),1);
					d.offset = d1.offset.*d2.rot + d2.offset;
				else																							% d2.rot is a matrix
					d.rot = d1.rot*d2.rot;
					d.offset = d1.offset*d2.rot + d2.offset;
				end
			end
			d.lablist_in = d1.lablist_in;
			w = prmapping(mfilename,'trained',d,getlabels(w2),getsize_in(w1),getsize_out(w2));
		else
			% Store a sequential mapping.
			w = sequential(w1,w2);
		end
		
	else  

		% Execution of the affine mapping.
		% R is a dataset or double, OFFSET defines the mapping.
		
		v = offset;
		[m,k] = size(R);
		d = +v;
		
		if all(size(v) == 0)
			d.rot = repmat(d.rot,1,k);
			d.offset = zeros(1,k);
		end

		if (size(d.rot,1) == 1) & (k > 1) 
			% No rotation, just a scaling
			x = zeros(m,k);
			Rdat = +R;
			if (m > k)							% Necessary switch for handling large feature sizes.
				for j=1:k
					x(:,j) = Rdat(:,j)*d.rot(j);
				end
			else
				for i=1:m
					x(i,:) = Rdat(i,:).*d.rot;
				end
			end
			x = x + repmat(d.offset,m,1);

		else		% Rotation.
			
			x = [+R,ones(m,1)] * [d.rot;d.offset];
			
		end
		
		if size(v,2) == 2 & size(x,2) == 1
			x = [x -x];
    end
		
    if isdataset(R)
      w = setdat(R,x,v);
    else
      w = x;
    end
		
	end

return;
