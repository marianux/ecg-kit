%BAGC Bag classifier for classifying sets of object instances
%
%		[WBAG,WOBJ] = BAGC(A,OBJCLASSF,BAGINDEX,BAGCOMBC,BAGCLASSF,BAGLAB)
%		 WBAG       = A*BAGC([],OBJCLASSF,BAGINDEX,BAGCOMBC,BAGCLASSF,BAGLAB)
%		 WBAG       = A*BAGC(OBJCLASSF,BAGINDEX,BAGCOMBC,BAGCLASSF,BAGLAB)
%		 D          = B*WBAG
%
% INPUT
%   A          Training dataset with object labels and bag indices 
%   B          Test Dataset with index list of bags, stored as label list 
%   OBJCLASSF  Trained or untrained object classifier, default QDC
%   BAGINDEX   String or a label_list_index defining in which label list
%              of A the bag indices are stored, default 2.
%   BAGCOMBC   Combiner for objects in a bag, default VOTEC
%   BAGCLASSF  Untrained classifier for bags, default FISHERC
%   BAGLAB     String or a label_list_index defining in which label list
%              of A the bag labels are stored. Objects with the same 
%              bag index should have the same bag label. 
%              Default is the current labeling of A
%
% OUTPUT
%   WBAG       Trained bag classifier
%   WOBJ       Trained object classifier
%   D          Classification matrix of bags in B
%
% DESCRIPTION
% This routine offers a classifier for bags (e.g. images) of objects
% (e.g. pixels) stored in a single dataset. The objects in the training 
% set A should have at least two labels: bag labels (the class of their bag)
% and bag indices, defining which objects belong to the same bag. These two
% label sets should be stored by the ADDLABELS command in the dataset A.
% Refer to the multi-labeling system (see MULTI_LABELING) offered by 
% PRTools. The current object labels of A can be the bag labels, but may
% also be different, e.g. true object labels.
%
% BAGINDEX should be a label_list_name or a label_list_index defining the 
% label list used for storing the bag indices that refer to the bag an 
% object belongs to. The same label_list_name or label_list_index should be
% used for defining the bags of the test objects in B.
%
% All objects in A are used to train the object classifier OBJCLASSF if it 
% is untrained. The current object labels are used for that. Classification 
% results of the objects in the same bag are combined by BAGCOMBC, which 
% can be any of the fixed combiners MEANC, PRODC, PERC, MAXC, etcetera.
% This results for every bag in a single confidence vector for the classes. 
%
% If an untrained bag classifier BAGCLASSF is supplied, the bag confidence
% vectors are used to train a bag classifier.
%
% New bags, organised in a dataset like B, with the proper bag indices per
% object stored in a label list with the same name or label_list_index as 
% used in A, can be classified by the bag classifier WBAG. 
%
% If no bag classifier BAGCLASSF was defined during training, just the
% results of the object classifier WOBJ are returned combined by BAGCOMBC
% over the objects in the same bag in B. In this case the final result is
% identical to B*(A*WOBJ)*BAGCC([],BAGCOMBC), provided that A has class
% labels and B is labeled by its bag indices.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, MULTI_LABELING, BAGCC, LOSO,
% DATASET/ADDLABELS, DATASET/CHANGELABLIST

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [out1,out2] = bagc(varargin)

	argin = shiftargin(varargin,'prmapping');
  argin = setdefaults(argin,[],qdc,2,votec,[],[]);
  if mapping_task(argin,'definition')
    out1 = define_mapping(argin,'untrained');
    out1 = setname(out1,'Set classifier');
  else
    [a,objclassf,bagindex,bagcombc,bagclassf,baglab] = deal(argin{:});
    if isuntrained(objclassf) | nargin > 2 | ~strcmp(getmapping_file(objclassf),mfilename)
      % train the mapping (classifier)

      % we need datasets with at least 1 object per class and 2 classes
      isvaldset(a,1,2);

      % if the object classifier is untrained, train it, else use it
      if isuntrained(objclassf)
        wobj = a*objclassf;
      else
        wobj = objclassf;
      end

      if ismapping(bagclassf) & isuntrained(bagclassf)

        % if the bag labels are not given, 
        % use the objcts labels for the bags too
        if isempty(setlab), setlab = curlablist(a); end

        % classifiy the dataset and change labeling to bag index
        x = changelablist(a*wobj,setindex);

        % avoid empty bags
        x = setlablist(x);

        % combine object results to bag results
        d = bagcc(x,bagcombc);

        % change to bag labels
        d = changelablist(d,baglab);

        % train bag classifier
        bagclassf = d*bagclassf;

        % get outputlabels
        labels_out = getlabels(bagclassf);

      else
        labels_out = getlabels(wobj);
      end

      % store all what is needed for execution in the mapping
      out1 = prmapping(mfilename,'trained',{wobj,bagcombc,bagclassf,bagindex,baglab}, ...
        labels_out, size(a,2),size(labels_out,1));

      % prevent batch execution
      out1 = setbatch(out1,0);

      % return the object classifier as well
      out2 = wobj;

    else % here we are for execution

      % the mapping is stored in objclassf
      w = getdata(objclassf);

      % save current lablist
      curlist = curlablist(a);

      % use the bag index for the test set if supplied
      if ~isempty(w{4}), testset = changelablist(a,w{4}); end

      % classify test set by the object classifier
      d = testset*w{1};

      % avoid empty bags
      d = setlablist(d);

      % combine objects in the same bag
      d = bagcc(d,w{2}); 

      % reset lablist for classification matrix
      d = changelablist(d,curlist);

      % apply the set classifier, if defined
      if ~isempty(w{3}), d = d*w{3}; end

      % that is it, define class labels as feature labels
      out1 = setfeatlab(d,getlabels(objclassf));

    end
    
  end
	
	
	
	
	
	