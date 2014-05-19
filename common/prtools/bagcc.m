%BAGCC Combining classifier for classifying bags of objects 
%
%		DBAG = BAGCC(DOBJ,COMBC)
%		DBAG = DOBJ*BAGCC([],COMBC)
%
% INPUT
%   DOBJ   Dataset, classification matrix, output of some base classifier
%   COMBC  Combiner, e.g. MAXC (default VOTEC)
%
% OUTPUT
%   DBAG   Dataset, classification matrix for the bags in DOBJ
%
% DESCRIPTION
% This routine combines object classification results of bags of objects 
% stored in DOBJ. It is assumed that the current labels of DOBJ are bag 
% identifiers and defining objects belonging to the same bag. Objects of 
% the same bag are combined by COMBC into a single classification result 
% and returned by DBAG. 
%
% DBAG gets as many objects as there are bags defined for DOBJ. Effectively
% the first object of every bag in DOBJ is replaced by the combined result
% and other objects of that bag are deleted. DBAG has the same feature 
% labels (most likely the class names) as DOBJ and stores as object 
% identifiers the bag identifiers stored in the label list of DOBJ. 
% A possible multi-labeling definition of DOBJ is preserved.
%
% This routine is called by BAGC where needed.
%
% SEE ALSO
% DATASETS, BAGC, MULTI-LABELING

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [dbag,id] = bagcc(dobj,combc)

if nargin < 2 | isempty(combc), combc = votec; end

if nargin < 1 | isempty(dobj)
	% define the mapping
	dbag = prmapping(mfilename,'untrained',{combc});
	dbag = setname(dbag,'Bag combiner');
else % execution
	
	% we should have a proper dataset
	isdataset(dobj);
	
	% the class names are the bag indentifiers
	bagnames = classnames(dobj);
	
	% retrieve datasize, and number of sets c
	[m,k,c] = getsize(dobj);
	
	% get number of objects for every set
	s = classsizes(dobj);
	
	% dobj is a classification matrix, so its features point to classes
	featlab = getfeatlab(dobj);
	
	% reserve spave for the result
	dbag = prdataset(zeros(c,k));

	% space the object identifiers of the first object per bag
	id = zeros(c,1);
	
	t = sprintf('Combining %i bags: ',c);
	prwaitbar(c,t);
	
	% run over all bags
	for j=1:c
		prwaitbar(c,j,[t int2str(j)]);
		
		% get the objects in the bag
		y = seldat(dobj,j);
		
		% the identifier of the first object
		id(j) = getident(y(1,:));
		
		%create a dataset with all objects in the bag concatenated horizontally
		y = +y';
		y = prdataset(y(:)');
		
		% give the the proper class labels
		y = setfeatlab(y,repmat(featlab,s(j),1));
		
		% now classifier combiners can be used
		dbag(j,:) = y*combc;
		
	end
	prwaitbar(0);
	
	% find the first objects of every set
	J = findident(dobj,id); 
	
	% and replace them by the bag combining result
	% so object labels become bag labels
	dbag = setdata(dobj(J,:),dbag);
	
	% give columns the classnames
	dbag = setfeatlab(dbag,featlab);
	
	% use set bag names as bag identifiers
	dbag = setident(dbag,bagnames);
end