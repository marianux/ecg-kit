%WEAKC Weak Classifier
%
%   [W,V] = WEAKC(A,ALF,ITER,CLASSF)
%   [W,V] = A*WEAKC([],ALF,ITER,CLASSF)
%   [W,V] = A*WEAKC(ALF,ITER,CLASSF)
%
% INPUT
%   A       Dataset
%   ALF     Fraction or number of objects to be used for training, see
%           GENDAT. Default: one object per class. For ALF is integer, ALF
%           objects per class are generated.
%   ITER    Number of trials
%   CLASSF  untrained classifier, default NMC
%
% OUTPUT
%   W       Best classifier over ITER runs
%   V       Cell array of all classifiers
%           Use VC = stacked(V) for combining
%   VC      Combined set of classifiers
%
% DESCRIPTION
% WEAKC uses subsampled versions of A for training. Testing is done
% on the entire training set A. The best classifier is returned in W.
%
%    VC = WEAKC(A,ALF,ITER,CLASSF,1)
%
% Combines all classifiers as a stacked combiner in VC.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, NMC, GENDAT, STACKED

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,v] = weakc(varargin)
	
%               INITIALISATION

	argin = shiftargin(varargin,'scalar',1);
  argin = setdefaults(argin,[],1,1,0,0);
  
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained','Weak');
    
%                 TRAINING

  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [a,n,iter,r,s] = deal(argin{:});

    if isscalar(n) && n >= 1
      n = n*ones(1,getsize(a,3));
    end
    v = {};
    emin = 1;

    for it = 1:iter              % Loop
      b = gendat(a,n);           % subsample training set
      if ~ismapping(r)           % select classifier and train
        if r == 0                % be consistent with old classfier selection
          ww = nmc(b); 
        elseif r == 1
          ww = fisherc(b); 
        elseif r == 2
          ww = udc(b);
        elseif r == 3
          ww = qdc(b);
        else
          error('Illegal classifier requested')
        end
      else
        if ~isuntrained(r)
          error('Input classifier should be untrained')
        end
        ww = b*r;
      end
      v = {v{:} ww};              % store all classifiers in v
                                  % select best classifier and store in w
      e = a*ww*testc;
      if e < emin
        emin = e;
        w = ww;
      end
    end

    if s == 1
      w = stacked(v);
    end
    
  end

return