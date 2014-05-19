%PRGLOBAL Define or reset PRTools globals
%
%   PRGLOBAL
%   PRGLOBAL RESET
%   PRGLOBAL(PAR1,VAL1,PAR2,VAL2, ....)
%
% INPUT
%   PAR1,PAR2  String with variable name of global variable
%   VAL1,VAL2  Desired value of global variable
%
% DESCRIPTION
% The functioning of PRTools depends on the values of a number of global
% variables. Many of them have their own routine to change them. This
% routine is a central routine to list them (no arguments), reset them
% to their initial values or set them to some user preference.
%
% The following global variables can be changed:
% PRMEMORY,  see corresponding routine
% DEFAULTBATCHSIZE,  see SETBATCH
% GRIDSIZE,  see corresponding routine and PLOTC
% CHECK_SIZES,  TRUE or FALSE, disables or enables some size checking
% STAMP_MAP,  TRUE or FALSE, see corresponding routine
% PRWAITBAR,  ON or OFF, see corresponding routine
% PRWARNING,  set level, see corresponding routine
% REGOPT_NFOLDS,  see REGOPTC
% REGOPT_ITERMAX,  see REGOPTC
% REGOPT_REPS,  see REGOPTC
%
% For more information just type PRGLOBAL
%
% SEE ALSO
% PRMEMORY, SETBATCH, GRIDSIZE, STAMP_MAP, PRWAITBAR, PRWARNING, REGOPTC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function prglobal(varargin)

global CHECK_SIZES
global REGOPT_NFOLDS REGOPT_REPS REGOPT_ITERMAX REGOPT_OPTCRIT 
global DEFAULTBATCHSIZE
if isempty(DEFAULTBATCHSIZE)
  DEFAULTBATCHSIZE = 1000;
end
 
if nargin == 0
  s = '\n          ';
  fprintf(['\n%5i  PRMEMORY/1e6, the maximum size of matrices intermediately created by' ...
    s 'PRTools. If set too high Matlab may crash. If set low PRTools may' s ...
    'generate many loops. \n'],prmemory/1e6);
  fprintf(['\n%5i  DEFAULTBATCHSIZE, the batch size set as default by setbatch. If set ' ...
    s 'too high Matlab may crash. If set low PRTools may generate many loops. \n'], ...
    DEFAULTBATCHSIZE);
  fprintf(['\n%5i  GRIDSIZE, to be used in scatterplots. Its value determines the' ...
    s 'accuracy by which classifiers are drawn using plotc. The default is' s ...
    '30 but sometimes 300 or more is needed causing long computing times.\n'],gridsize);
  fprintf(['\n%5i  CHECK_SIZES, true (1) or false(0), determines whether in sequential' ...
    s 'combining of datassets and mappings output and input sizes should be' s ...
    'checked on matching.\n'],CHECK_SIZES);
  fprintf(['\n%5i  STAMP_MAP, disabled (0), retrieval only (1), storage and retrieval (2)' ...
    s 'of trained mappings to avoid recomputation of precomputed results.\n'],stamp_map);
  fprintf(['\n%5s  PRWAITBAR, the status of prwaitbar, ''on'' or ''off''\n'],prwaitbar);
  fprintf(['\n%5i  PRWARNING, the level of prwarning.\n'],prwarning);
  fprintf(['\n%5i  REGOPT_NFOLDS, the number of folds used in the crossvalidation' ...
    s 'procedure by which the classifier performance is estimated for a' s ...
    'single setting of the parameters during their optimization by REGOPTC.\n'],REGOPT_NFOLDS);
  fprintf(['\n%5i  REGOPT_ITERMAX, the maximum number of iterations in the optimization' ...
    s 'procedure by which the classifier performance is estimated for a' s ...
    'single setting of the parameters during their optimization by REGOPTC.\n'],REGOPT_ITERMAX);
  fprintf(['\n%5i  REGOPT_REPS, the number of repetitions used in the crossvalidation' ...
    s 'procedure by which the classifier performance is estimated for a' s ...
    'single setting of the parameters during their optimization by REGOPTC.\n'],REGOPT_REPS);
  fprintf(['\n          Note that automatic parameter optimization for classifiers may slow' ...
    s 'down training by a factor REGOPT_NFOLDS*REGOPT_ITERMAX*REGOPT_REPS.\n']);
  
elseif (nargin == 1) & strcmpi(varargin{1},'reset')
  prmemory(50000000);
  gridsize(30);
  prwarning(1);
  prwaitbar;
  CHECK_SIZES = true;
  stamp_map(0);
  REGOPT_NFOLDS = 5;
  REGOPT_REPS = 1;
  REGOPT_ITERMAX = 20;
  DEFAULTBATCHSIZE = 1000;
elseif nargin ~= 2*floor(nargin/2)
  error('Wrong number of input arguments')
else
  for j=1:2:nargin
    switch lower(varargin{j})
      case 'globalprmemory'
        prmemory(varargin{j+1});
      case 'defaultbatchsize'
        DEFAULTBATCHSIZE = varargin{j+1};
      case 'current_gridsize'
        gridsize(varargin{j+1});
      case 'prwaitbar'
        prwaitbar(varargin{j+1});
      case 'prwarning'
        prwarning(varargin{j+1});
      case 'check_sizes'
        if all(varargin{j+1} ~= [true,false])
          error('CHECK_SIZES should either be ''true'' or ''false''')
        end
        CHECK_SIZES = varargin{j+1};
      case 'stamp_map'
        stamp_map(varargin{j+1});
      case 'regopt_nfolds'
        REGOPT_NFOLDS = varargin{j+1};
      case 'regopt_reps'
        REGOPT_REPS = varargin{j+1};
      case 'regopt_itermax'
        REGOPT_ITERMAX = varargin{j+1};
      otherwise
        error('Unknown global parameter')
    end
  end
end
        
        
    