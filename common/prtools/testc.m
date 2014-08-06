%TESTC Test classifier, error / performance estimation
%
%   [E,C] = TESTC(A*W,TYPE)
%   [E,C] = TESTC(A,W,TYPE)
%    E    = A*W*TESTC([],TYPE)
%
%   [E,F] = TESTC(A*W,TYPE,LABEL)
%   [E,F] = TESTC(A,W,TYPE,LABEL)
%    E    = A*W*TESTC([],TYPE,LABEL)
%
% INPUT
%   A      Dataset
%   W      Trained classifier mapping
%   TYPE   Type of performance estimate, default: probability of error
%   LABEL  Target class, default: none
% 
% OUTPUT
%   E      Error / performance estimate
%   C      Number of erroneously classified objects per class.
%          They are sorted according to A.LABLIST
%   F      Error / performance estimate of the non-target classes
% 
% DESCRIPTION
% This routine supplies several performance estimates for a trained
% classifier W based on a test dataset A. It is possible to supply a cell  
% array of datasets {A*W}, or a cell array of datasets {A} or a cell array
% of classifiers {W}. In case A as well as W is a cell array, W might be
% 2-dimensional in with as many columns as A has datasets. See DISPERROR
% for an example.
%
% A should contain test objects for every class assigned by W.
% Objects in A belonging to different classes than defined for W as well
% as unlabeled objects are neglected. Note that this implies that TESTC
% applied to a rejecting classifier (e.g. REJECTC) estimates the
% performance on the not rejected objects only. By 
%   [E,C] = TESTC(A,W); E = (C./CLASSSIZES(A))*GETPRIOR(A)';
% the classification error with respect to all objects in A may be
% computed. Use CONFMAT for an overview of the total class assignment 
% including the unlabeled (rejected) objects.
%
% In case of missing classes in A, [E,C] = TESTC(A*W) returns in E a NaN
% but in C still the number of erroneously classified objects per class.
%
% If LABEL is given, the performance estimate relates just to that class as
% target class. If LABEL is not given a class average is returned weighted
% by the class priors.
%
% The following performance measures are supported for TYPE:
% 'crisp'       Expected classification error based on error counting,
%               weighted by the class priors (default).
% 'FN'          E False negative
%               F False positive
% 'TP'          E True positive
%               F True negative
% 'soft'        Expected classification error based on soft error
%               summation, i.e. a sum of the absolute difference between
%               classifier output and target, weighted by class priors.
% 'F'           Lissack and Fu error estimate
% 'mse'         Expected mean square difference between classifier output
%               and target (based on soft labels), weighted by class
%               priors.
% 'auc'         Area under the ROC curve (this is an error and not a
%               performance!). For multi class problems this is the
%               weigthed average (by class priors) of the one-against-rest
%               contributions of the classes.
% 'precision'   E Fraction of true target objects among the objects
%                 classified as target. The target class is defined by LABEL. 
%                 Priors are not used.
%               F Recall, fraction of correctly classified objects in the 
%                 target class. Priors are not used.
% 'sensitivity' E Fraction of correctly classified objects in the target
%                 class (defined by LABEL). Priors are not used. 
%                 Sensitivity as used her is identical to recall.
%               F Specificity, fraction non target objects that are not 
%                 classified into the target class (defined by LABEL). 
%                 Priors are not used.
%
% EXAMPLES
% See PREX_PLOTC.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, CONFMAT, REJECTC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: testc.m,v 1.19 2010/02/18 15:57:07 duin Exp $

function [OUT1,OUT2] = testc(a,w,type,label)

  if nargin < 4, label = []; end
	if nargin < 3, type = [];  end
	if nargin < 2, w = []; end
	if nargin < 1, a = []; end
	
	if isstr(w) % takes care of testc(a*w,type,label)
		label = type; type = w; w = [];
	end
% 		
% 	if isempty(type)
% 		type = 'crisp';
% 	end
	
	if (isempty(a)) % prepares a*testc([],w,type,label), or a*testc([],type,label)
		out1 = prmapping(mfilename,'fixed_cell',{w,type,label});
		out1 = setname(out1,'testc');
		out1 = setbatch(out1,0); % Don't run in batch mode
		out2 = [];
		
	elseif (~ismapping(w) & ~iscell(w)) | (ismapping(w) & isfixed(w) & strcmp(getname(w),'testc'))
		% call like testc(a*w,type,label), or a*testc([],w,type,label) which
		% results in testc(a*w,V) in which V = testc([],type,label)

		if ismapping(w) % retrieve parameters stored in testc([],w,type,label
			label = getdata(w,3);
			type = getdata(w,2);
			w = getdata(w,1);
		end
	
		if (iscell(a))
  		% If this argument is a cell array, recursively call this
  		% function to get errors for all elements in the cell array.
			out1 = zeros(size(a));
			out2 = cell(size(a));
     % if isempty(w)
        
			for j1 = 1:size(a,1)
				for j2 = 1:size(a,2)
					[out1(j1,j2),out2{j1,j2}] = feval(mfilename,a{j1,j2},w,type,label);
				end
			end
			
    elseif (isdatafile(a)) % datafile needs some handling as we need to 
			                     % process all objects separately
	    c = getsize(a,3);
	    out2 = zeros(1,c);
	    next = 1;
      a = setprior(a,getprior(a));
	    while next > 0
		    [b,next] = readdatafile(a,next);
				if isempty(w)
					[out1,class_err] = feval(mfilename,prdataset(b));
				else
		    	[out1,class_err] = feval(mfilename,b,w,type,label);
				end
		    out2 = out2 + class_err;
	    end
	    if isempty(a.prior)
		    out1 = sum(out2)/size(a,1);
	    else
		    p = getprior(a);  
				csizes = classsizes(a);
				if any(csizes == 0)
					out1 = NaN;
					prwarning(1,'Some classses have no test objects')
				else
		    	out1 = (out2./classsizes(a))*p';
				end
	    end
      
		else % here we are for the real work
			
			isdataset(a); % we need a real dataset for evaluation
			if (ismapping(w) & istrained(w))
				a = a*w;
			end
			
			fflab = renumlab(getfeatlab(a));
			if any(fflab == 0) % reject option! Allowed?
				if isempty(strmatch(type,char('crisp','FN','TP','precision','sensitivity')))
					error('Reject option incompatible with desired error measure')
				else % remove all objects to be rejected
					reject_col = find(fflab == 0);
					[ma,J] = max(+a,[],2);
					L = find(J~=reject_col);
					a = a(L,:);
				end
			end
			a = a*maxc; % takes care that every class appears as a single column in a
			%a = remclass(a);
			lablist = getlablist(a);  % classes of interest
			featlab = getfeatlab(a);
			flab = renumlab(featlab,lablist);

			csizes = classsizes(a);
 			if any(flab == 0) 
				prwarning(1,'Some classes assigned by the classifier have no test objects')
			end
			if any (csizes == 0)
				if nargout < 2 | ~isempty(label) % no error / performance measure can be returned
					error('Some classes have no test objects')
				else % we can, however, return, the error per class
					prwarning(2,'Some classes have no test objects')
	        c = getsize(a,3);
          I = matchlablist(a*labeld,lablist);
          nlab = getnlab(a);
					OUT2 = zeros(1,c);
					for j=1:c
            J = find(nlab==j);
						OUT2(j) = sum(I(J)~=j);
					end
					OUT1 = NaN;
				end
				return
			end
			
			clab = renumlab(lablist,featlab);
			if any(clab == 0)
				prwarning(1,'Superfluous test objects found, they will be neglected')
				J = find(clab~=0);
				a = seldat(a,J);
				a = setlablist(a);
        csizes = csizes(J);
			end
			
			[m,k,c] = getsize(a); p = getprior(a); labtype = getlabtype(a);
      
			if isempty(type) % set default error measure types
				if islabtype(a,'crisp')
					type = 'crisp';
				elseif islabtype(a,'soft')
					type = 'soft';
				else
					type = 'mse'
				end
			end

			confm = cmat(a); % compute confusion matrix
			if isempty(label)
				lablist = getlablist(a);
				out2 = (csizes - diag(confm)');
				out = zeros(1,c);
				out1 = 0;
				for j = 1:c			% compute crit one-against-rest	
% 					confm2 = [confm(j,j) sum(confm(j,:))-confm(j,j); sum(confm(:,j))-confm(j,j) ...
% 						sum(confm(:))-sum(confm(:,j))-sum(confm(j,:)) + confm(j,j)];
					confm2 = [confm(j,j) csizes(j)-confm(j,j); sum(confm(:,j))-confm(j,j) ...
						sum(confm(:))-sum(confm(:,j))-csizes(j) + confm(j,j)];
					b = seldat(a,j);
					out(j) = comp_crit(type,confm2,a,j,lablist(j,:));
					if isempty(a.prior) 
						out1 = out1 + out(j);
					else
						out1 = out1 + p(j) * out(j) / size(b,1);
					end
				end
				if isempty(a.prior)
					out1 = out1 / m;
				end
			else
				n = getclassi(a,label);
				confm2 = [confm(n,n) sum(confm(n,:))-confm(n,n); sum(confm(:,n))-confm(n,n) ...
						sum(confm(:))-sum(confm(:,n))-sum(confm(n,:)) + confm(n,n)];
				[out1,out2] = comp_crit(type,confm2,a,n,label);
				out1 = out1/csizes(n);
				out2 = out2/(m-csizes(n));
			end
			
		end

	elseif (iscell(a)) || (iscell(w))
		
		% If there are two input arguments and either of them is a cell array,
		% call this function on each of the cells.

		% Non-cell array inputs are turned into 1 x 1 cell arrays.

		if (~iscell(a)), a = {a}; end
		if (~iscell(w)), w = {w}; end
    
    if size(a,2) > 1 & size(w,1) == 1
      % repeat w for all input datasets
      w = repmat(w,size(a,2),1);
    end
    
		if (min(size(a) > 1))
			error('2D cell arrays of datasets not supported')
    end

		% Now call this function for each combination of 
		% dataset A{I} and mappings W{I,J}.

		out1 = cell(numel(a),size(w,2)); 
		out2 = out1;
		for i=1:numel(a)
			for j=1:size(w,2)
				[out1{i,j},out2{i,j}] = feval(mfilename,a{i}*w{i,j},type,label);
			end
		end

	else

		% Assert that the second argument is a trained mapping, and call
		% this function on the mapped data.
		ismapping(w); istrained(w);
		[out1,out2]= feval(mfilename,a*w,type,label);

	end

	% If there are no output arguments, display the error(s) calculated.
	% Otherwise, copy the calculated errors to the output arguments.

	if (nargout == 0) & (nargin > 0)
		if (iscell(a))
			if (nargin == 1) || (isempty(w) && isempty(type) && isempty(label))
				for j1 = 1:size(a,1)
					for j2 = 1:size(a,2)
            fprintf(1,'%6.4f %s on %15s\n',out1(j1,j2),...
              getname(a{j1,j2}),getuser(a{j1,j2},'evaluated_by'));
						%disp(['Mean classification error on ' ...
            %num2str(size(a{j1,j2},1)) ' test objects: ' num2str(out1(j1,j2))]);
					end
				end
			else
				fprintf('\n  Test results for');
				disperror(a,w(1,:),cell2mat(out1));
			end
		else
			if (nargin == 1)
				disp(['Mean classification error on ' num2str(size(a,1)) ' test objects: ' num2str(out1)])
			else
				if ~isempty(w) %DXD empty mapping can happen after a*w*testc
					fprintf(' %s',getname(w,20));
				else  %DXD is this a good alternative?
					fprintf(' %s',getname(a,20));
				end
				fprintf(' %5.3f',out1);
				fprintf('\n');
			end
		end
	else
		OUT1 = out1;
		OUT2 = out2;
	end

return

%TESTAUC Multiclass error area under the ROC
%
%   E = TESTAUC(A*W)
%   E = TESTAUC(A,W)
%   E = A*W*TESTAUC
%
% INPUT
%   A  Dataset to be classified
%   W  Classifier
%
% OUTPUT
%   E  Error, Area under the ROC
%
% DESCRIPTION
% The area under the error ROC is computed for the datset A w.r.t. the
% classifer W. The estimator is based on a rank analysis of the classifier
% outcomes. Ties are broken by a two-way sorting and averaging. 
%
% The multiclass situation is solved by averaging over all outcomes of
% the one-against-rest ROCs.
%
% Note that E is an error and not a performance measure like the AUC often
% used in literature.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, TESTC, PRROC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function e = testauc(a,w,n)

	  if nargin < 3, n = []; end
  
	if (nargin == 0) | (isempty(a))
		% No input arguments given: return mapping information.
		e = prmapping(mfilename,'fixed',{label});
		e = setbatch(e,0);
    return
	elseif (nargin == 1 | isempty(w))
		% Classification matrix already computed
    d = a;
  else
    % Compute classification matrix now
    d = a*w;
  end

  [m,k,c] = getsize(d);
  s = classsizes(d);

  if k == 1 % classifier with a single class outcome, make two for consistency
    d = [d 1-d];
    k = 2;
  end
  
  if isempty(n)
    e = zeros(1,c);
    for j = 1:c
      e(j) = auc_one(d,s,j);
    end
    e = e*getprior(d)';
  else
    e = auc_one(d,s,n);
  end
  
return

function e = auc_one(d,s,j)

% compute AUC for class j versus rest
 
  m = size(d,1);
  lablist = getlablist(d);   % class names
  n = findfeatlab(d,lablist(j,:));
  % forward sorting
  [ds,J1] = sort(-d(:,n));  [j1,K1] = sort(J1);
  % backward sorting to solve ties
  [ds,J2] = sort(flipud(-d(:,n))); [j2,K2] = sort(J2); K2 = flipud(K2);
  % get all object indices for this class
  K = findnlab(d,j);
  % retrieve number of wrong pairs
  e = (sum(K1(K)) + sum(K2(K))-(s(j)*(s(j)+1)))/2;
  % error contribution
  e = e / ((m-s(j))*s(j));
  
return

function confm = cmat(a)

% simplified confusion matrix procedure, class order as in c.lablist
% a should be a classification matrix with the same feature labels
% (no doubles) as a.lablist

lablist = getlablist(a);
featlab = getfeatlab(a);
N = getsize(a,3);
flab = renumlab(featlab,lablist);
nlab = getnlab(a);
aa = +a;

confm = zeros(N,N);
for j=1:N
	J = find(nlab==j);
	[mx,K] = max(aa(J,:),[],2);
	confm(j,:) = histc(flab(K)',1:N);
end

function [out1,out2] = comp_crit(type,c,a,n,label)

% c : 2 x 2 confusion matrix
% a : classification data 
% n : relevant class

switch type
	case 'crisp'
		out1 = c(1,2);
		out2 = c(2,1);
	case 'FN'
		out1 = c(1,2);
		out2 = c(2,1); % FP
	case 'TP'
		out1 = c(1,1);
		out2 = c(2,2); % TN
	case 'precision'
		out1 = c(1,1)/(c(1,1)+c(2,1));
		out1 = out1*(c(1,1)+c(1,2)); % sum of per sample contributions
		out2 = c(1,1)/(c(1,1)+c(1,2)); % recall (=sensitivity)
		out2 = out2*(c(2,2)+c(2,1)); % sum of per sample contributions
	case 'sensitivity'
		out1 = c(1,1)/(c(1,1)+c(1,2)); 
		out1 = out1*(c(1,1)+c(1,2)); % sum of per sample contributions
		out2 = c(2,2)/(c(2,1)+c(2,2)); % specificity
		out2 = out2*(c(2,2)+c(2,1)); % sum of per sample contributions
	case 'soft'  % normalised difference between desired and real targets
		a = setlabtype(a,'soft')*classc;
		t = gettargets(a);
		k = findfeatlab(a,label);
		d = abs(+a(:,k) - t(:,n));
		J = find(isnan(d));
		d(J) = ones(size(J));
		out1 = sum(d)/2; % needed for consistency as every error is counted twice
		%out1 = sum(d);
		out2 = [];
	case 'F'     % Lissack and Fu error
		b = seldat(a,n)*classc;
		out1 = sum(1-max(+b,[],2));
		out2 = [];
	case {'mse','MSE'}
		k = findfeatlab(a,label);
		b = seldat(a,n);
		out1 = sum((+b(:,k)-gettargets(b)).^2);
	case {'nmse','NMSE'} % use normalised outputs
		k = findfeatlab(a,label);
		b = seldat(a,n)*classc;
		out1 = sum((+b(:,k)-gettargets(b)).^2);
	case {'auc','AUC'}
		out1 = testauc(a,[],n)*(c(1,1)+c(1,2));% sum of per sample contributions
		out2 = [];
	otherwise
		error('Error / performance type not found')
end
		
		
