%SETBATCH Enables/diasables execution in batch mode of mappings
%
%    W = SETBATCH(W,FLAG,BATCHSIZE,OBJSIZE)
%    W = W*SETBATCH([],FLAG,BATCHSIZE,OBJSIZE)
%
% INPUT
%    W          Mapping
%    FLAG       Switches batch processing on (TRUE, default) or off (FALSE)
%    BATCHSIZE  Number of objects per batch (default 1000)
%    OBJSIZE    Number of objects above which batch processing is enabled
%               (default BATCHSIZE)
%
% OUTPUT
%    W          Mapping
%
%DESCRIPTION
%Mappings might be processed in batch mode in case of large datasets or
%datafiles. This may solve problems in computing A*W (A is a datase/file) 
%by executing the result in batches if intermediate memory demand is too
%high. In case the final dataset is too large this will not be of use as it
%should fit in memory anyway. 
%
%For some mappings batch processing is not possible as the processing of
%objects is not fully independent. For that reason it is skipped during
%training of mappings and classifiers.
%
%SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%MAPPINGS, GETBATCH

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function w = setbatch(varargin)

  global DEFAULTBATCHSIZE
  if isempty(DEFAULTBATCHSIZE)
    DEFAULTBATCHSIZE = 1000;
  end
  argin = setdefaults(varargin,[],true,DEFAULTBATCHSIZE,[]);
  [w,n,batchsize,objsize] = deal(argin{:});
  if isempty(objsize), objsize = batchsize; end
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'combiner');
    w = setuser(w,{n,batchsize,objsize},'batch');
  else
    if isempty(objsize), objsize = batchsize; end
    if all([0,1] ~= n)
      error('Flag should be 0 or 1')
    end
    ismapping(w);
    w = setuser(w,{n,batchsize,objsize},'batch');
  end
  
return
