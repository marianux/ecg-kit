%PROXM Proximity mapping
% 
%  W = PROXM(A,TYPE,P,WEIGHTS)
%  W = A*PROXM([],TYPE,P,WEIGHTS)
%  W = A*PROXM(TYPE,P,WEIGHTS)
%  D = B*W
% 
% INPUT
%  A       Dataset used for representation
%  B       Dataset applied to the representation set
%  TYPE    Type of the proximity (optional; default: 'distance')
%  P       Parameter of the proximity (optional; default: 1)
%  WEIGHTS Weights (optional; default: all 1)
%
% OUTPUT
%  W       Proximity mapping
%  D       Dataset, proximity matrix between B and A
% 
% DESCRIPTION  
% Computation of the [K x M] proximity mapping (or kernel) defined by 
% the [M x K] dataset A. Unlabeled objects in A are neglected.
% If B is an [N x K] dataset, then D=B*W is the [N x M] proximity matrix
% between B and A. The proximities can be defined by the following types
% (A and B are here 1 x K vectors): 
% 
% 	'polynomial'   | 'p': SIGN(A*B'+1).*(A*B'+1).^P
% 	'homogeneous'  | 'h': SIGN(A*B').*(A*B').^P
% 	'exponential'  | 'e': EXP(-(||A-B||)/P)
% 	'radial_basis' | 'r': EXP(-(||A-B||.^2)/(P*P))
% 	'sigmoid'      | 's': SIGM((SIGN(A*B').*(A*B'))/P)
% 	'distance'     | 'd': ||A-B||.^P
% 	'minkowski'    | 'm': SUM(|A-B|.^P).^(1/P)
% 	'city-block'   | 'c': SUM(|A-B|)
%   'cosine'       | 'o': 1 - (A*B')/||A||*||B||
%   'ndiff'        | 'n': SUM(A~=B)
%   'hellinger'    | 'g': ||A.^0.5-B.^0.5||
% 
% In the polynomial case for a non-integer P, the proximity is computed 
% as D = SIGN(S+1).*ABS(S+1).^P, in order to avoid problems with negative
% inner products S = A*B'. The features of the objects in A and B may be 
% weighted by the weights in the vector WEIGHTS.
%
% Note that for computing (squared) Euclidean distances DISTM is much
% faster than PROXM.
%
% EXAMPLES
% W = A*proxm('m',1);              % define L1 distances
% W = A*proxm('m',1); D = B*W;     % L1 distances between B and A
% W = proxm('r',2)*mapex; D = A*W; % Distances between A and itself
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, DISTM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% 27-11-2008, Marc Hulsman, MR1
% Added kernel normalization.

% $Id: proxm.m,v 1.10 2010/01/15 13:39:01 duin Exp $

function W = proxm(varargin)
		 
  varargin = shiftargin(varargin,'char');
  argin = setdefaults(varargin,[],'d',1,[]);
  if mapping_task(argin,'definition')
    W = define_mapping(argin,'untrained','Proximity mapping');
    return
  end
  
  [A,type,s,weights] = deal(argin{:});
	[m,k] = size(A);
  
	if (isstr(type))
		% Definition of the mapping, just store parameters.

		all = char('polynomial','p','homogeneous','h','exponential','e','radial_basis','r', ...
	  		       'sigmoid','s','distance','d','minkowski','m','city-block','c',...
                   'cosine','o','uniform','u','ndiff','n','hellinger','g');
		if (~any(strcmp(cellstr(all),type)))
			error(['Unknown proximity type: ' type])
		end
		A = cdats(A,1);
		[m,k] = size(A);
		%if (isa(A,'prdataset'))
		if (isdataset(A))
			A = testdatasize(A);
			W = prmapping(mfilename,'trained',{A,type,s,weights},getlab(A), ...
                            getfeatsize(A),getobjsize(A));
		else
			W = prmapping(mfilename,'trained',{A,type,s,weights},[],k,m);
		end
		W = setname(W,'Proximity mapping');
										   
	elseif isa(type,'prmapping')  
		
		% Execution of the mapping: compute a proximity matrix
		% between the input data A and B; output stored in W.
		w = type;
		[B,type,s,weights] = deal(w.data{1},w.data{2},w.data{3},w.data{4});
		[kk,n] = size(w);
		if (k ~= kk)
			error('Mismatch in sizes between the mapping and its data stored internally.'); 
		end
		if (~isempty(weights))		
			if (length(weights) ~= k), 
				error('Weight vector has a wrong length.'); 
			end
			A = A.*(ones(m,1)*weights(:)');
			B = B.*(ones(n,1)*weights(:)');
    end

    testdatasize(m*n);
		switch type
			case {'polynomial','p'}
				D = +(A*B');
				D = D + ones(m,n);
				if (s ~= round(s))
					D = +sign(D).*abs(D).^s;
				elseif (s ~= 1)
					D = D.^s;
				end
	
			case {'homogeneous','h'}
				D = +(A*B'); 
				if (s ~= round(s))
					D = +sign(D).*abs(D).^s;
				elseif (s ~= 1)
					D = D.^s;
				end
	
			case {'sigmoid','s'}
				D = +(A*B'); 
				D = sigm(D/s);
		
			case {'city-block','c'}
				D = zeros(m,n);
				t = sprintf('Processing %i objects: ',n);
				prwaitbar(n,t);
				for j=1:n
					prwaitbar(n,j,[t int2str(j)]);
					D(:,j) = sum(abs(A - repmat(+B(j,:),m,1)),2);
				end
				prwaitbar(0);
		
			case {'minkowski','m'}
				D = zeros(m,n);
				t = sprintf('Minkowski distance for %i objects: ',n);
				prwaitbar(n,t);
				for j=1:n
					prwaitbar(n,j,[t int2str(j)]);
					D(:,j) = sum(abs(A - repmat(+B(j,:),m,1)).^s,2).^(1/s);
				end
				prwaitbar(0);
		
			case {'exponential','e'}
				D = dist2(B,A);
				D = exp(-sqrt(D)/s);
		
			case {'radial_basis','r'}
				D = dist2(B,A);
				D = exp(-D/(s*s));
		
			case {'distance','d'}
				D = dist2(B,A);
				if s ~= 2
					D = sqrt(D).^s;
				end
                
			case {'cosine','o'}
				D = +(A*B'); 
				lA = sqrt(sum(A.*A,2));
				lB = sqrt(sum(B.*B,2));
				D = 1 - D./(lA*lB');
				
			case {'ndiff','n'}
				D = zeros(m,n);
				t = sprintf('Processing %i objects: ',n);
				prwaitbar(n,t);
				for j=1:n
					prwaitbar(n,j,[t int2str(j)]);
					D(:,j) = sum(+A ~= repmat(+B(j,:),m,1),2);
				end
				prwaitbar(0);
		
			case {'hellinger','g'}
				D = dist2(sqrt(B),sqrt(A));
        D = sqrt(D);
				
			otherwise
				error('Unknown proximity type')
    end
    if isdataset(A)
      W = setdat(A,D,w);	
    else
      W = D;
    end
		
	else
		error('Illegal TYPE argument.')
	end

return;



function D = dist2(B,A)
% Computes square Euclidean distance, avoiding large matrices for high 
% dimensional data

	ma = size(A,1);
	mb = size(B,1);
	
	if isdatafile(A)
	
		D = zeros(ma,mb);
		next = 1;
		while next > 0
			[a,next,J] = readdatafile(A,next);
			D(J,:) = +dist2(a,B);
		end
		
	elseif isdatafile(B)

		D = zeros(ma,mb);
		next = 1;
		while next > 0  % we need another version of readdatafile here, as due
			[b,next,J] = readdatafile2(B,next); % to persistent variables double
			D(:,J) = +dist2(A,b); % looping can not be handled correctly
		end
		
	else % A and B are not datafiles

		D = ones(ma,1)*sum(B'.*B',1); 
		D = D + sum(A.*A,2)*ones(1,mb);
		D = D - 2 .* (+A)*(+B)';
		% Clean result.
		J = find(D<0);
		D(J) = zeros(size(J));
		
	end
	
return
