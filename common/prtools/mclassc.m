%MCLASSC Computation of multi-class classifier from 2-class discriminants
%
%  W = MCLASSC(A,CLASSF,MODE)
%  W = A*MCLASSC([],CLASSF,MODE)
%  W = A*MCLASSC(CLASSF,MODE)
%
% INPUT
%   A       Dataset
%   CLASSF  Untrained classifier
%   MODE    Type of handling multi-class problems (optional; default: 'single')
%
% OUTPUT
%   W       Combined classifier
%
% DESCRIPTION
% For default MODE = 'single', the untrained classifier CLASSF is called to
% compute C classifiers between each of the C classes in the dataset A and
% the remaining C-1 classes. The result is stored into the combined
% classifier W.
%
% For MODE = 'multi' the untrained classifier CLASSF is trained between all
% pairs of classes as well as between each class and all other classes.
% This total set of C classes is combined by MINC. The use of soft labels
% is supported.
%
% EXAMPLES
% W = MCLASSC(GENDATM(100),QDC,'MULTI');
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, MINC.

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: mclassc.m,v 1.9 2009/08/18 23:09:01 duin Exp $

function [w,varargout] = mclassc(varargin)

  varargout = repmat({[]},[1, max((nargout-1),0)]);
  argin = shiftargin(varargin,'prmapping');
  argin = setdefaults(argin,[],[],'single');
  
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained');
    return
  end
    
  [a,classf,mode] = deal(argin{:});	
	if ~isa(classf,'prmapping') || ~isuntrained(classf)
		error('Second parameter should be untrained mapping')
	end

	islabtype(a,'crisp','soft');
	isvaldfile(a,1,2); % at least 1 object per class, 2 classes

	[m,k,c] = getsize(a);
	
	if c == 2
		[w,varargout] = prmap(a,classf);
		return
	end

	varout = {};
  lablist = getlablist(a);
	s = sprintf('Multi-class: %4i classifiers: ',c); 
	prwaitbar(c,s);

	switch mode

	 case 'single'
	  w = [];
	  %	lablist = getlablist(a);
		
	  for i=1:c
			prwaitbar(c,i,[s int2str(i)]);
		  if islabtype(a,'crisp')
			  mlab = 2 - (getnlab(a) == i);
			  aa = setlabels(a,mlab);
				aa = remclass(aa);  % remove empty classes
				%aa = setnlab(a,mlab);
				%aa = setlablist(aa,[1 2]');
				if ~isempty(a.prior)
					aa = setprior(aa,[a.prior(i),1-a.prior(i)]');
				end
		  elseif islabtype(a,'soft')
			  atargets = gettargets(a);
			  targets = [atargets(:,i) 1-atargets(:,i)]; %assumes soft labels sum to one???
			  aa = prdataset(+a,mlab,targets,'lablist',[1 2]');
		  end
			varo = varargout;
      [v,varo{:}] = prmap(aa,classf);
			varout = [varout; varo];
			w = [w,setlabels(v(:,1),lablist(i,:))];
	  end

	 case 'multi'
	  w = [];

		nclassf = 0;
	  nlab = getnlab(a);
	  for i1=1:c
			prwaitbar(c,i1,[s int2str(i1)]);
		  lab = lablist(i1,:);
		  
		  J1 = find(nlab==i1);
		  if islabtype(a,'crisp')
			  mlab = ones(m,1);
			  mlab(J1) = zeros(length(J1),1);
			  aa = setlabels(a,mlab);
				aa = remclass(aa); % remove empty classes
		  else
			  problab = gettargets(a);
			  mlab = [problab(:,i1) sum(problab,2)-problab(:,i1)];
			  aa = settargets(a,mlab,[1 2]');
		  end		
		  I1 = [1:c]; I1(i1) = [];
			varo = varargout;
      [v,varo{:}] = prmap(aa,classf);
			varout = [varout; varo];
      w = [w,setlabels(v(:,1),lab)];

		  for i2 = I1 
			  if islabtype(a,'crisp')
				  J2 = find(nlab==i2);
				  v = aa([J1;J2],:)*classf;
			  else
				  mlab2 = problab(:,[i1 i2]);
				  v = setlabels(aa,mlab2)*classf;
			  end
			  w = [w,setlabels(v(:,1),lab)];
				nclassf = nclassf+1;
		  end
	  end
	  w = minc(w);

	 otherwise
	  error('Unknown mode')
	end
	prwaitbar(0);	
	
	w = setname(w,getname(classf));
	w = setsize(w,[k,c]);
	w = setcost(w,a);

	if ~isempty(varout)
    varargout = num2cell(varout',2)';
  end
  
	return
