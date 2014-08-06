%WVOTEC Weighted combiner (Adaboost weights)
%
%  W = WVOTEC(A,V)   compute weigths and store
%  W = WVOTEC(V,U)   Construct weighted combiner using weights U
%
%  INPUT
%    A      Labeled dataset
%    V      Parallel or stacked set of trained classifiers
%    U      Set of classifier weights
%
%  OUTPUT
%    W      Combined classifier
%
% DESCRIPTION
% The set of trained classifiers V is combined using weighted
% majority voting. If given the weights U are used. If not
% given, the weights are computed from the classification
% results of the labeled dataset A using 0.5*log((1-E)/E)
% if E is the classifier error. 
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS,

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = wvotec(a,v)

if nargin < 1 | isempty(a)
	%w = prmapping(mfilename,'untrained');
	w = prmapping(mfilename,'combiner');
elseif nargin < 2
	if ismapping(a) 
		w = prmapping(mfilename,'untrained',{a}); % a is the set of base classifiers
	else
		error('Illegal call')
	end
elseif isdataset(a) | isdouble(a)      % train or classify
	if isuntrained(v)
		v = a*v;								           % train base classifiers
	end
	if ~strcmp(v.mapping_file,mfilename) % training (base classifiers are already trained)
		if isparallel(v)                   % parallel combiner
			n = 0;
			e = zeros(1,length(v.data));
			for j=1:length(v.data)
				vv = v.data{j};
				d = a(:,n+1:n+size(vv,1))*vv*classc;
				e(j) = testmc(d);
				n = n+size(vv,1);
			end
		elseif isstacked(v)                % stacked combiner
			e = zeros(1,length(v.data));
			for j=1:length(v.data)
				vv = v.data{j};
				%e(j) = testc(a,vv,'soft');
				e(j) = testc(a,vv);
			end
		else
			error('Classifier combination expected')
		end
		                  % Find weights								
		L = find(e < 1-max(getprior(a))); % take classifier better than prior
		alf = zeros(1,length(e));
		alf(L) = log((1-e(L))./e(L))/2;
		alf = alf/sum(alf);
		                  % construct the classifier
		[m,k,c] = getsize(a);
		w = prmapping(mfilename,'trained',{v,alf},getlabels(vv),k,c);
		w = setname(w,'Weighted Voting');
	else                                 % testing
		alf = v.data{2};                   % get the weights
		u = v.data{1};                     % get the set of classifiers
		m = size(a,1);
		dtot = zeros(m,size(v,2));
		if isparallel(u)                   % parallel combiner
			n = 0;
			for j=1:length(u.data)           % weight them 
				vv = u.data{j};
				aa = a(:,n+1:n+size(vv,1));
				d = a(:,n+1:n+size(vv,1))*vv;
				[dd,jj] = max(+d,[],2);
				dd = zeros(size(dtot));
				dd([1:m]'+(jj-1)*m) = alf(j);
				dtot = dtot + dd;
				n = n+size(vv,1);
			end
		elseif isstacked(u)                % stacked combiner
			for j=1:length(u.data)           % weight them
				vv = u.data{j};
				d = a*vv;
				[dd,jj] = max(+d,[],2);
				dd = zeros(size(dtot));
				dd([1:m]'+(jj-1)*m) = alf(j);
				dtot = dtot + dd;
			end
		else
			error('Classifier combination expected')
		end
		w = setdat(d,dtot);
  end
  
else                  % store classifier from given weights
    
  ismapping(a);
  u = v;              % the weights
  v = a;              % the combined classifier
  n = length(v.data);
  if length(u) ~= n
    error('Wrong number of weights given')
  end
  [k,c] = getsize(v.data{1});
	w = prmapping(mfilename,'trained',{v,u},getlabels(v{1}),k,c);
	w = setname(w,'Weighted Voting');
end

		
			