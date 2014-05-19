%KERNELM Kernel mapping, dissimilarity representation
% 
%   [W,J] = KERNELM(A,KERNEL,SELECT,P1,P2 , ...)
%    W    = A*KERNELM([],KERNEL,SELECT,P1,P2 , ...)
%    K    = B*W
%
% INPUT
%   A,B         Datasets
%   KERNEL      Untrained kernel / dissimilarity representation,
%               a mapping computing proximities between objects.
%               default: Euclidean dissimilarities: PROXM([],'d',1)
%   SELECT      Name of object selection procedure, see below
%   P1,P2, ...  Additional parameters for SELECT
%
% OUTPUT
%   W          Mapping
%   J          Vector with indices of selected objects for representation
%   K          Kernel matrix, dissimilarity representation, 
%              size [SIZE(B,1) LENGTH(J)]
%
% DESCRIPTION
% Computes the kernel mapping W for the representation objects in A. The 
% computation of the kernel matrix, which is a proximity matrix (similarities
% or dissimilarities) should be defined in KERNEL by an untrained mapping like
% PROXM for predefined proximities or USERKERNEL for user specified
% proximities.
% A*KERNEL should 'train' the kernel, i.e. specify A as representation set.
% B*(A*KERNEL) should compute the kernel matrix: a dataset.
%
% Initially, the kernel mapping has a size [SIZE(A,2) SIZE(A,1)]. For
% increased efficiency or accuracy the representation set may be reduced
% by a routine given by the string SELECT to select to objects J, using 
% possibly additional parameters P1, P2, etcetera. This option of
% representation set reduction is the only difference between the use of
% KERNELM and routines like PROXM and USERKERNEL.
%
% The following choices for SELECT are supported:
% 
% 'random'    random selection of P1 objects, maximum P2
% 'gendat'    [X,Y,J] = GENDAT(A,P1)
% 'kcentres'  [LAB,J] = KCENTRES(DISTM(A),P1,P2)
% 'modeseek'  [LAB,J] = MODESEEK(DISTM(A),P1)
% 'edicon'    J = EDICON(DISTM(A),P1,P2,P3)
% 'featsel'   J = +FEATSELM(A*KERNELM(A,TYPE,P),P1,P2,P3)
%
% REFERENCES
% 1. E.Pekalska, R.P.W.Duin, P.Paclik, Prototype selection for dissimilarity-
% based classification, Pattern Recognition, vol. 39, no. 2, 2006, 189-208.
% 2. E.Pekalska and R.P.W.Duin, The Dissimilarity Representation for Pattern
% Recognition, Foundations and Applications, World Scientific, 2005, 1-607.
% 
% EXAMPLE
% A = GENDATB;
% W = (SCALEM*KERNELM([],[],'random',5)*LOGLC); 
% SCATTERD(A)
% PLOTC(A*W)
%
% SEE ALSO
% DATASETS, MAPPINGS, PROXM, USERKERNEL

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: kernelm.m,v 1.7 2007/07/10 08:25:29 duin Exp $

function w = kernelm(a,kernel,select,varargin);

		
	if isempty(varargin), varargin = {[] [] []}; end
	if length(varargin) == 1, varargin = {varargin{:} [] []}; end
	if length(varargin) == 2, varargin = {varargin{:} []}; end
	
	if (nargin < 3) 
		select = []; 
		prwarning(4,'No representation set reduction specified.')
	end
	if (nargin < 2) | isempty(kernel)
		prwarning(3,'No kernel mapping specified. Euclidean distances assumed')
		kernel = proxm([],'d',1); 
	end

	if (nargin < 1) | (isempty(a))  
		% Definition: an untrained mapping.
		w = prmapping(mfilename,{kernel,select,varargin{:}});
		w = setname(w,'Kernel mapping');
		
	elseif isstr(kernel) % old format of call: kernelm(a,type,p,n), training
		type = kernel;
		p = select;
		if ~isempty(varargin)
			if length(varargin) > 1
				error('Wrong parameters supplied')
			end
			n = varargin{1};
		else
			n = [];
		end
		[m,k] = size(a);
		kernel = proxm([],type,p);
		w = prmapping(mfilename,'trained',{a*kernel},getlab(a),k,m);
		if ~isempty(n)
			w = w*pcam(a*w,n);
		end
		w = setname(w,'Kernel Mapping');

	elseif isa(kernel,'prmapping') & ~strcmp(getmapping_file(kernel),mfilename) % training
	
		a = testdatasize(a);
		a = testdatasize(a,'objects');
		isuntrained(kernel);
		[m,k] = size(a);
		%w = prmapping('kernelm','trained',{a*kernel},getlab(a),k,m);
		if isempty(select)
			w = prmapping('kernelm','trained',{a*kernel},getlab(a),k,m);
		elseif ismapping(select)
			r = a*select;
			w = prmapping('kernelm','trained',{r*kernel},getlab(a),k,size(r,1));
		else
			switch select
				case 'random'
					J = randperm(m);
					n = varargin{1};
					if isempty(n) | n > m
						error('Number of objects to be selected not given or too large')
					end
					if n < 1, n = ceil(n*m); end % fraction given
					if ~isempty(varargin{2})
						n = min(n,varargin{2});
					end
					J = J(1:n);
				case 'gendat'    
					[x,y,J] = gendat(a,varargin{1});
				case 'kcentres'  
					[lab,J] = kcentres(distm(a),varargin{1:2});
				case 'modeseek'  
					[lab,J] = modeseek(distm(a),varargin{1});
				case 'edicon'    
					J = edicon(distm(a),varargin{1:3});
				case 'featsel'
					w = prmapping('kernelm','trained',{a*kernel},getlab(a),k,m);
					J = +featselm(a*w,varargin{1:3});
				otherwise
					error('Unknown choice for object selection')
			end
			% redefine mapping with reduced representation set
			labels_out = getlab(a);
			w = prmapping('kernelm','trained',{a(J,:)*kernel},labels_out(J,:),k,length(J));
		end
		w = setname(w,'Kernel Mapping');
		
	else % Execution of the mapping, w will be a dataset.
		
		kern = getdata(kernel,1); % trained kernel is stored in datafield
		K = a*kern;
		w = setdat(a,K,kern);
		
	end

return;
