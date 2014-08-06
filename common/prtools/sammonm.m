%SAMMONM Sammon mapping
%
%   W = SAMMONM(A,K,MAX)
%   W = A*SAMMONM([],K,MAX)
%   W = A*SAMMONM(K,MAX)
%   D = B*W
%
% INPUT
%   A    Dataset, used for training the mapping
%   B    Dataset, same dimensionality as A, to be mapped
%   K    Target dimension of mapping (default 2)
%   MAX  Maximum number of iterations, default 100
%
% OUTPUT
%   W    Trained mapping
%   D    K-dimensional dataset
%
% DESCRIPTION
% This is a simplified interface to the more complex MDS routine for high
% dimensional data visualisation. The output is a non-linear projection of
% the original vector space to a K-dimensional target space. 
%
% The main differences with MDS are that SAMMONM operates on feature based
% datasets, while MDS expects dissimilarity matrices; MDS maps new objects
% by a second optimization procedures minimizing the stress for the test
% objects, while SAMMONM uses a linear mapping between dissimilarities and
% the target space. See also PREX_MDS for examples. A different procedure
% for the same purpose is TSNEM.
%
% EXAMPLE
% prdatasets;            % make sure prdatasets is in the path
% a = satellite;         % 36D dataset, 6 classes, 6435 objects
% [x,y] = gendat(a,0.5); % split in train and test set
% w = x*sammonm;         % compute mapping
% figure; scattern(x*w); % show trainset mapped to 2D: somewhat overtrained
% figure; scattern((x+randn(size(x))*1e-5)*w): % some noise helps
% figure; scattern(y*w); % show test set mapped to 2D
%
% REFERENCES
% 1. JW Sammon Jr A nonlinear mapping for data structure analysis,
%   IEEE Transactions on Computers C-18, pp. 401-409,1969.
% 2. E. Pekalska, D. de Ridder, R.P.W. Duin, and M.A. Kraaijveld, A new
%    method of generalizing Sammon mapping with application to algorithm
%    speed-up, ASCI99, Proc. 5th Annual ASCI Conf., 1999, 221-228. [<a href="http://rduin.nl/papers/asci_99_mds.pdf">pdf</a>]
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PCAM, MDS, TSNEM, PREX_MDS, SCATTERD, SCATTERN

% Copyright: E. Pekalska, R.P.W. Duin, r.p.w.duin@37steps.com

function out = sammonm(varargin)

  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],2,100,0.01);
  if mapping_task(argin,'definition')
    out = define_mapping(argin,'untrained');
  elseif mapping_task(argin,'training')
    [a,k,max_iter,pow] = deal(argin{:});
    d = sqrt(distm(+a));
    opt.maxiter = max_iter;
    w = mds(d,a*pcam(a,k),opt);
    y = getdata(w,1);
    % refer to distances with a low power to get some smoothness
    v = prpinv(d.^pow)*y;
    out = trained_mapping(a,{v,a,pow},2);
  elseif mapping_task(argin,'trained execution')
    [b,w] = deal(argin{1:2});
    [v,a,pow] = getdata(w);
    out = setdata(b,(sqrt(distm(b,a)).^pow)*v);
  else
    error('Illegal call');
  end

return