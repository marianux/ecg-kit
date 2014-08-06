%GAUSSM Trainable mapping, mixture of Gaussians (MoG) density estimate
%
%  W = GAUSSM(A,K,R,S,M)
%  W = A*GAUSSM([],K,R,S,M)
%  W = A*GAUSSM(K,R,S,M)
%
% INPUT
%   A     Dataset
%   K     Number of Gaussians per class
%   R,S,M Regularization parameters, 0 <= R,S <= 1, see QDC
%
% OUTPUT
%   W     Mixture of Gaussians density estimate
%
% DESCRIPTION
% Estimation of a PDF by the dataset A by a mixture of Gaussians procedure.
% Use is made of EMCLUST(A,QDC,K). Unlabeled objects are neglected, unless
% A is entirely unlabeled or double. Then all objects are used. If A is a
% multi-class crisp labeled dataset the densities are estimated class by
% class and then weighted and combined according their prior probabilities.
% Use +A instead of A to obtain a single set of Gaussians. In all cases,
% just single density estimator W is returned.
%
% Note that it is necessary to set the label type of A to soft labels
% (A = LABTYPE(A,'soft') in order to use the traditional EM algorithm
% based on posterior probabilities instead of using crisp labels.
% 
% The mapping W may be applied to a new dataset B using DENSITY = B*W.
%
%  W = A*GAUSSM
%
% uses a single Gaussian per class (K=1) and no regularisation. If
% regulariisation is desired, also K should be supplied.
%
% EXAMPLE
% a = gendatb;
% w = a*gaussm(2);
% scatterd(a)
% plotm(w)
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, QDC, MOGC, EMCLUST, PLOTM, TESTC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: gaussm.m,v 1.8 2009/02/02 21:57:39 duin Exp $

function w = gaussm(varargin)

  mapname = 'Mixture of Gaussians';
	argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],1,0,0,[]);
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained');
    w = setname(w,mapname);
  else
    [a,n,r,s,dim] = deal(argin{:});
	
    if isa(a,'prdataset')
      labname = getname(a);
    else
      labname = '';
    end
    
    if ((~isdataset(a) && ~isdatafile(a)) || ...
        (getsize(a,3) ~= 1 && islabtype(a,'crisp')))
      % this handles the multiclass situation
      w = mclassm(a,prmapping(mfilename,n),'weight');
      w = setlabels(w,labname);
      w = setname(w,mapname);
      
    else
      % here we have a single class
      [m,k] = getsize(a);
      if n == 1
        % just one Gaussian, esitmate it
        [U,G] = meancov(a);
        res.mean = +U;
        res.cov  = G;
        res.prior= 1;
        w = normal_map(res,labname,k,1);
      else
        % multiple Gaussians, run EM
        [e,v] = emclust(a,qdc([],r,s,dim),n);	
        ncomp0 = size(v.data.mean,1);	
        iter = 0;
        while (ncomp0 ~= n & iter < 5)   % repeat until exactly n components are found
          [e,v1] = emclust(a,qdc([],r,s,m),n);				
          ncomp1 = size(v1.data.mean,1);
          if ncomp1 > ncomp0
            v = v1;
            ncomp0 = ncomp1;
          end
          iter = iter + 1;
        end
        res = v.data;
        res.nlab = ones(n,1); % defines that all Gaussian components have to be
                              % combined into a single class.
        w = prmapping('normal_map','trained',res,labname,k,1);
      end
      w = setname(w,mapname);
      
    end
    
  end

return
