%FIXEDCC Construction of fixed combiners, back-end routine
%
%   V = FIXEDCC(A,W,TYPE,NAME,PAR)
%
% INPUT
%   A      Dataset
%   W      A set of classifier mappings
%   TYPE   String defining the type of combination rule
%   NAME   The name of this combination rule, arbitrary
%   PAR    Possible parameter for combiner defined by TYPE
%
% OUTPUT
%   V      Mapping
%
% DESCRIPTION
% Define a mapping V which applies the combination rule TYPE to the
% set of mappings W. The set of mappings W should be a parallel
% combination (see MAPPINGS).
%
% TYPE defines the combining rule and can be any of the following:
% average, min, max, mean, median, prod, vote,
%
% Note that average is only possible for affine 2-class classifiers.
%
% When W is a set of classifiers and A a dataset (possibly the result
% of B*W, where W is again a set of classifiers) then:
%
%   V = FIXEDCC(W,[],TYPE)   combines the mappings W with the comb. rule TYPE
%   V = FIXEDCC(A,[],TYPE)   computes the combining output for dataset A
%   V = FIXEDCC(A,W,TYPE)    computes the combining output for dataset A,
%                            where W is trained using A
% 
% This is a back-end routine. Users should directly call the fixed combiners.
%
% EXAMPLES
% See prex_combining.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, VOTEC, MAXC, MEANC, MEDIANC, MINC, PRODC, AVERAGEC

% $Id: fixedcc.m,v 1.3 2009/02/04 10:53:10 duin Exp $

function v = fixedcc(a,w,type,name,par)

		if nargin == 0
		% Without input, just return an completely empty combiner mapping (even
		% without defining the combination rule):

		v = prmapping(mfilename,'combiner');
		return
	end	

	if nargin < 5, par = []; end
	if isempty(a)
		% just return a combiner mapping with just the combining rule
		
		v = prmapping(mfilename,'combiner',{[],type,name,par});
		
	elseif isa(a,'prmapping') & isempty(w)
		
		% v = comb_classifier*fixedcc,
		% we combine a set of classifiers, not interesting, except for the
		% 'average' type, there a new affine transformation is computed:
		
		if isuntrained(a) % store untrained comb_classifier
			
			v = prmapping(mfilename,'untrained',{a,type,name,par});
			
		else              % handle or store trained comb_classifier and all info
			
			[nclass,classlist] = renumlab(getlabels(a));

			% Special case is to average affine coefficients for 'average' combining
			switch type
			 case 'average'
			  if ~(isaffine(a) & size(classlist,1) == 2)
				  error('Average combining only possible for affine 2-class classifiers')
			  end
			  n = length(a.data);
			  rot = zeros(size(a,1),1);
			  off = 0;
			  for j=1:length(a.data)
				  rot = rot + a.data{j}.data.rot(:,1);
				  off = off + a.data{j}.data.offset(1);
			  end
			  v = affine(rot/n,off/n);			
			otherwise
			  % Standard procedure: make a new trained mapping
			  v = prmapping(mfilename,'trained',{a,type,name,par});
			end
			v = set(v,'size_in',size(a,1),'size_out',size(classlist,1),'labels',classlist);
			v = setcost(v,a);
		end	
		
	elseif isa(a,'prdataset') & isempty(w) 
		
		% Call like v = dataset*fixedcc,
		% Here the work will be done:
		% the dataset has already been mapped through the mappings, and the outputs
		% should be processed according to the combiner type.
		
		% get all the relevant parameters:
		[m,k] = size(a);
		featlist = getfeatlab(a);
		if isempty(featlist)
			prwarning(2,'No class names given: numbering inserted')
			nclass = [1:k]';
			classlist = [1:k]';
		else
			[nclass,classlist] = renumlab(featlist);
		end
		c = size(classlist,1);
		d = zeros(size(a,1),c);
		b = +a;  % the classifier outputs to be processed
		%DXD: I need for my one-class classifiers that the feature domains
		%are retained:
		newfeatdom = getfeatdom(a);
		if ~isempty(newfeatdom)
			newfeatdom = newfeatdom(1:size(d,2));
		end

		% for each of the classes the outputs should now be combined to a new
		% one, using the combining rule:
		for j=1:c
			
			J = find(nclass==j);
			
			switch type
				
			case 'min'
			  d(:,j) = min(b(:,J),[],2);
			  
			case 'max'
			  d(:,j) = max(b(:,J),[],2);
			  
			case 'mean'
			  d(:,j) = mean(b(:,J),2);
			  
			case 'prod'
			  %d(:,j) = prod(b(:,J),2);
				d(:,j) = exp(sum(log(b(:,J)),2));
				
			case 'median'
			  d(:,j) = median(b(:,J),2);
			  
			case 'vote' % Assumes that classifier outcomes are well ordered in b
				% For voting we cannot combine the basic classifier outputs,
				% but we should use the classifier labels:
				n = size(a,2) / c;
				if ~isint(n)
					error('All classifiers should refer to all classes')
				end
				% First get the votes for each of the classes:
				[dummy,fl] = renumlab(featlist);
				mlab = zeros(m,n);
        % we are in a loop over all clases, but we treat here all classes
        % simultaneously and break below. First run over all classifiers
				for j=1:n
					J = [(j-1)*c+1:j*c];
					labels = labeld(a(:,J));
					[dummy,nlab,ll] = renumlab(fl,labels);
					mlab(:,j) = nlab;
				end
				% Then count the number of votes for every class
				for j=1:c
					d(:,j) = (sum(mlab==j,2)+1)/(n+c);
				end
				%DXD: for this voting rule, the feature domain will change:
				for k=1:c
					newfeatdom{k} = [0 inf; 0 inf];
        end
        % and we are done
				break
        
			case 'perc'
				d(:,j) = prctile(b(:,J),par,2);
				
			case 'average'
				error([newline 'Average combiner should directly call the classifiers' ...
				newline 'e.g. A*AVERAGEC([W1 W2 W3]), or A*([W1 W2 W3]*AVERAGEC)'])
				     
			otherwise
			  error(['Unknown method for fixed combining: ' type])
			end
			
		end
		
		v = setdata(a,d,classlist);
		%DXD: I need for my one-class classifiers that the feature domains
		%are retained:
		if ~isempty(newfeatdom)
			v = setfeatdom(v,newfeatdom);
		end
		
		
	elseif (isa(a,'double') || isa(a,'prdataset')) && isa(w,'prmapping')
		
		% call like v = dataset * trained combiner (e.g. a*votec([u v w]))

		% This means that we first have to map the data through the mappings, and
		% then map this new dataset through the combiner:

		if strcmp(getmapping_file(w),mfilename)
			% Then we already have a nice combining rule:
			% get the relevant parameters:
			type = w.data{2};
			name = w.data{3};
			par - w.data{4};
			% Evaluate the mapped data (a*w.data{1}) by this combining rule:
			v = feval(mfilename,a*w.data{1},[],type,name,par);
		else
			% We will use the parameters given in the argument list
			% evaluate the mapped data (a*w) by this combining rule:
			v = feval(mfilename,a*w,[],type,name,par);
		end
		
	else                   % this should not happen
		
		error('Call cannot be parsed')
		
	end

	if isa(v,'prmapping')
		v = setname(v,name);
	end

	return
