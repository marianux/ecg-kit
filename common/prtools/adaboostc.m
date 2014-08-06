%ADABOOSTC
%
% [W,V,ALF] =  ADABOOSTC(A,CLASSF,N,RULE,VERBOSE);
%
% INPUT
%   A       Dataset
%   CLASSF  Untrained weak classifier
%   N       Number of classifiers to be trained
%   RULE    Combining rule (default: weighted voting)
%   VERBOSE Suppress progress report if 0 (default)
%
% OUTPUT
%   W       Combined trained classifier
%   V       Cell array of all classifiers
%           Use VC = stacked(V) for combining
%   ALF     Weights
%
% DESCRIPTION
%
% Computation of a combined classifier according to adaboost.
%
% In total N weighted versions of the training set A are generated
% iteratevely and used for the training of the specified classifier.
% Weights, to be used for the probabilities of the objects in the training
% set to be selected, are updated according to the Adaboost rule.
%
% The entire set of generated classifiers is given in V.
% The set of classifier weigths, according to Adaboost is returned in ALF
%
% Various aggregating possibilities can be given in 
% the final parameter rule:
% []:      WVOTEC, weighted voting.
% VOTEC    voting
% MEANC    sum rule
% AVERAGEC averaging of coeffients (for linear combiners)
% PRODC    product rule
% MAXC     maximum rule
% MINC     minimum rule
% MEDIANC  median rule
%
% REFERENCE
% Ji Zhu, Saharon Rosset, Hui Zhou and Trevor Hastie, 
% Multiclass Adaboost. A multiclass generalization of the Adaboost 
% algorithm, based on a generalization of the exponential loss.
% http://www-stat.stanford.edu/~hastie/Papers/samme.pdf
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
% (Multiclass correction by Marcin Budka, Bournemouth Univ., UK)

%function [W,V,alf] = adaboostc(a,clasf,n,rule,verbose)
function [out,V,alf] = adaboostc(varargin)

%% INITIALISATION
argin = setdefaults(varargin,[],nmc,100,[],0);
if mapping_task(argin,'definition')
  
  out = define_mapping(argin,'untrained','Adaboost');
  
%% TRAINING
elseif mapping_task(argin,'training')
  
  [a,clasf,n,rule,verbose] = deal(argin{:});
  [m,k,c] = getsize(a);
  V = [];
  laba = getlab(a);
  u = ones(m,1)/m;			% initialise object weights
  alf = zeros(1,n);			% space for classifier weights
  isseparable = 0;          % check if we can make 0 error
  if verbose && k == 2
    figure(verbose);
    scatterd(a);
  end

  %% generate n classifiers
  for i = 1:n
    b = gendatw(a,u,m);             % sample training set
    b = setprior(b,getprior(a));	% use original priors
    w = b*clasf;                    % train weak classifier
    ra = a*w;                       % test weak classifier

    if verbose && k == 2
      plotc(w,1); drawnow
    end
	
    labc = labeld(ra);
    diff = sum(labc~=laba,2)~=0;	% objects erroneously classified
    erra = sum((diff).*u);          % weighted error on original dataset

    if (erra==0)
        isseparable = 1;
        V = w;
        break;
    end
    if (erra < (1-1/c))        % if classifier better then random guessing...
      alf(i) = 0.5*(log((1-erra)/erra) + log(c-1));
      correct = find(diff==0); % find correctly classified objects
      wrong = find(diff==1);   % find incorrectly classified objects
      u(correct) = u(correct)*exp(-alf(i));	% give them the ...
      u(wrong) = u(wrong)*exp(alf(i));	  	% proper weights
      u = u./sum(u);                        % normalize weights
    else
      alf(i) = 0;
    end
	
    if verbose
      disp([erra alf(i) sum(alf)])
    end
    V = [V w];                       % store all classifiers

  end

  %% combine and return
  if isseparable
      W = V;
      W = setname(W,['Boosted ',getname(V)]);
  else
    if isempty(rule)
        W = wvotec(V,alf);             % default is weighted combiner
    else
        W = traincc(a,V,rule);         % otherwise, use user supplied combiner
    end
  end

  if verbose > 0 && k == 2
    plotc(W,'r',3)
    ee = a*W*testc;
    title(['Error: ', num2str(ee)]);
  end
  
  out = W;

else
  error('Illegal call')
end

return
