%MDS Trainale mapping for multidimensional scaling, a variant of Sammon mapping 
%
%   [W,J,STRESS] = MDS(DT,Y,OPTIONS)
%   [W,J,STRESS] = MDS(DT,K,OPTIONS)
%   [W,J,STRESS] = DT*MDS([],K,OPTIONS)
%   [W,J,STRESS] = DT*MDS(K,OPTIONS)
%              X = DS*W
%    
% 
% INPUT
%   DT      Square (M x M) dissimilarity matrix used for training
%   DS      N x M dissimilarity matrix between testset and trainset
%   Y       M x K matrix containing starting target configuration, or
%   K       Desired output dimensionality (default 2)
%   OPTIONS Various parameters of the  minimization procedure put into 
%           a structure consisting of the following fields: 'q', 'optim',
%           'init','etol','maxiter', 'isratio', 'itmap', 'st' and 'inspect'
%           (default:
%
%              OPTIONS.q       = 0          
%              OPTIONS.optim   = 'pn'      
%              OPTIONS.init    = 'cs'
%              OPTIONS.etol    = 1e-6 (the precise value depends on q)
%              OPTIONS.maxiter = 100 
%              OPTIONS.isratio = 0 
%              OPTIONS.itmap   = 'yes' 
%              OPTIONS.st      = 1 
%              OPTIONS.inspect = 2).
% 
% OUTPUT
%   W       Multidimensional scaling mapping
%   J       Index of points removed before optimization
%   STRESS  Vector of stress values
%   X       N x K dataset with output configuration
%
% DESCRIPTION  
% Finds a nonlinear MDS map (a variant of the Sammon map) of objects
% represented by a symmetric distance matrix D with zero diagonal, given
% either the dimensionality K or the initial configuration Y. This is done
% in an iterative manner by minimizing the Sammon stress between the given
% dissimilarities D (DT or DS) and the distances DY in the K-dimensional  
% target space:
%
%   E = 1/(sum_{i<j} D_{ij}^(Q+2)) sum_{i<j} (D_{ij} - DY_{ij})^2 * D_{ij}^Q
%
% If D(i,j) = 0 for any different points i and j, then one of them is
% superfluous. The indices of these points are returned in J.
%
% There is a simplified interface to MDS, called SAMMONM. The main
% differences are that SAMMONM operates on feature based datasets, while MDS 
% expects dissimilarity matrices; MDS maps new objects by a second
% optimization procedures minimizing the stress for the test objects, while
% SAMMONM uses a linear mapping between dissimilarities and the target
% space. See also PREX_MDS for examples.
%
% OPTIONS is an optional variable, using which the parameters for the mapping
% can be controlled. It may contain the following fields: 
%   
%   Q        Stress measure to use (see above): -2,-1,0,1 or 2. 
%   INIT     Initialisation method for Y: 'randp', 'randnp', 'maxv', 'cs' 
%            or 'kl'. See MDS_INIT for an explanation.
%   OPTIM    Minimization procedure to use: 'pn' for Pseudo-Newton or
%            'scg' for Scaled Conjugate Gradients.
%   ETOL     Tolerance of the minimization procedure. Usually, it should be 
%   MAXITER  in the order of 1e-6. If MAXITER is given (see below), then the 
%            optimization is stopped either when the error drops below ETOL or 
%            MAXITER iterations are reached.
%   ISRATIO  Indicates whether a ratio MDS should be performed (1) or not (0).
%            If ISRATIO is 1, then instead of fitting the dissimilarities 
%            D_{ij}, A*D_{ij} is fitted in the stress function. The value A 
%            is estimated analytically in each iteration.
%   ITMAP    Determines the way new points are mapped, either in an iterative 
%            manner ('yes') by minimizing the stress; or by a linear projection 
%            ('no').
%   ST       Status, determines whether after each iteration the stress should 
%   INSPECT  be printed on the screen (1) or not (0). When INSPECT > 0, 
%            ST = 1 and the mapping is onto 2D or larger, then the progress 
%            is plotted during the minimization every INSPECT iterations.
% 
% Important:
% - It is assumed that D either is or approximates a Euclidean distance 
%     matrix, i.e: 
%            D_{ij} = sqrt (sum_k(x_i - x_j)^2). 
% - Missing values can be handled; they should be marked by 'NaN' in D. 
%
% EXAMPLES:
% opt.optim = 'scg';
% opt.init  = 'cs'; 
% D  = sqrt(distm(a)); % Compute the Euclidean distance dataset of A
% w1 = mds(D,2,opt);   % An MDS map onto 2D initialized by Classical Scaling,
%                      % optimized by a Scaled Conjugate Gradients algorithm
% n  = size(D,1);
% y  = rand(n,2);
% w2 = mds(D,y,opt);   % An MDS map onto 2D initialized by random vectors
%
% z = rand(n,n);       % Set around 40% of the random distances to NaN, i.e. 
% z = (z+z')/2;        % not used in the MDS mapping
% z = find(z <= 0.6);
% D(z) = NaN;
% D(1:n+1:n^2) = 0;    % Set the diagonal to zero
% opt.optim = 'pn';
% opt.init  = 'randnp'; 
% opt.etol  = 1e-8;    % Should be high, as only some distances are used
% w3 = mds(D,2,opt);   % An MDS map onto 2D initialized by a random projection
%
% REFERENCES
% 1. M.F. Moler, A Scaled Conjugate Gradient Algorithm for Fast Supervised
%    Learning', Neural Networks, vol. 6, 525-533, 1993.
% 2. W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery,
%    Numerical Recipes in C, Cambridge University Press, Cambridge, 1992. 
% 3. I. Borg and P. Groenen, Modern Multidimensional Scaling, Springer
%    Verlag, Berlin, 1997. 
% 4. T.F. Cox and M.A.A. Cox, Multidimensional Scaling, Chapman and Hall, 
%    London, 1994.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PREX_MDS, MDS_CS, SAMMONM, TSNEM

%
% Copyright: E. Pekalska, R.P.W. Duin, e.pekalska@37steps.com
%

function [w,J,err,opt,y] = mds(varargin)

  mapname = 'MDS';
	argin = shiftargin(varargin,{'scalar'});
  argin = setdefaults(argin,[],2,[]);
  
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained',mapname);
    
  else % Train, evaluate or extend a mapping.
  
    [D,y,options] = deal(argin{:});
    opt = mds_setopt(options);

  % YREP contains representative objects in the projected MDS space, i.e. 
  % for which the mapping exists. YREP is empty for the original MDS, since 
  %  no projection is available yet.

    yrep = [];              

    if (isdataset(D) | isa(D,'double'))
      [m,mm] = size(D);

      % Convert D to double, but retain labels in LAB.
      if (isdataset(D)), lab = getlab(D); D = +D; else, lab = ones(m,1); end;

      if (ismapping(y))

        % The MDS mapping exists; this means that YREP has already been stored.
        pars = getdata(y); [k,c] = size(y);

        y     = [];        % Empty, should now be found.
        yrep  = pars{1};  % There exists an MDS map, hence YREP is stored.
        opt   = pars{2};  % Options used for the mapping.
        II    = pars{3};  % Index of non-repeated points in YREP.
        winit = pars{4};  % The Classical Scaling map, if INIT = 'cs'.
        v     = pars{5};  % Weights used for adding new points if ITMAP = 'no'.
        n = c;            % Number of dimensions of the projected space.

        % Initialization by 'cs' is not possible when there is no winit 
        % (i.e. the CS map) and new points should be added.
        if (strcmp(opt.init,'cs')) & (isempty(winit))
          prwarning(2,'OPTIONS.init = cs is not possible when adding points; using kl.');
          opt.init = 'kl';  
        end

        % If YREP is a scalar, we have an empty mapping.
        if (max(size(yrep)) == 1)
          y    = yrep;
          yrep = [];
        end

        % Check whether D is a matrix with the zero diagonal for the existing map.
        if (m == mm) & (length(intersect(find(D(:)<eps),1:m+1:(m*mm))) >= m)
          w = yrep;         % D is the same matrix as in the projection process; 
          return            % YREP is then the solution
        end

        if (length(pars) < 6) | (isempty(pars{6}))
          yinit = [];
        else
          yinit = pars{6};   % Possible initial configuration for points to 
                            % be added to an existing map
          if (size(yinit,1) ~= size(D,1))
            prwarning(2,'the size of the initial configuration does not match that of the dissimilarity matrix, using random initialization.')
            yinit =[];
          end
        end

      else

        % No MDS mapping available yet; perform the checks.

        if (~issym(D,1e-12))
          prwarning(2,'D is not a symmetric matrix; will average.'); 
          D = (D+D')/2;
        end

        % Check the number of zeros on the diagonal

        if (any(abs(diag(D)) > 0))
          error('D should have a zero diagonal'); 
        end
      end

    else    % D is neither a dataset nor a matrix of doubles
      error('D should be a dataset or a matrix of doubles.');
    end

    if (~isempty(y))

      % Y is the initial configuration or N, no MDS map exists yet;  
      % D is a square matrix.

      % Remove identical points, i.e. points for which D(i,j) = 0 for i ~= j.
      % I contains the indices of points left for the MDS mapping, J those
      % of the removed points and P those of the points left in I which were
      % identical to those removed.

      [I,J,P] = mds_reppoints(D);
      D = D(I,I);           
      [ni,nc] = size(D);

      % NANID is an extra field in the OPTIONS structure, containing the indices
      % of NaN values (= missing values) in distance matrix D.

      opt.nanid = find(isnan(D(:)) == 1); 

      % Initialise Y.

      [m2,n] = size(y);
      if (max(m2,n) == 1)    % Y is a scalar, hence really is a dimensionality N.
        n = y;
        [y,winit] = mds_init(D,n,opt.init);
      else
        if (mm ~= m2)
          error('The matrix D and the starting configuration Y should have the same number of columns.');
        end
        winit = [];
        y = +y(I,:);
      end

      % The final number of true distances is:
      no_dist = (ni*(ni-1) - length(opt.nanid))/2;

    else                      

      % This happens only if we add extra points to an existing MDS map. 
      % Remove identical points, i.e. points for which D(i,j) = 0 for i ~= j.
      % I contains the indices of points left for the MDS mapping, J those
      % of the removed points and P those of the points left in I which were
      % identical to those removed.

      [I,J,P] = mds_reppoints(D(:,II));
      D = D(I,II);             
      [ni,nc] = size(D);
      yrep = yrep(II,:);     
      n = size(yrep,2);     

      % NANID is an extra field in the OPTIONS structure, containing the indices
      % of NaN values (= missing values) in distance matrix D.

      opt.nanid = find(isnan(D(:))); 

      % Initialise Y. if the new points should be added in an iterative manner.

      [m2,n] = size(yrep);           
      if (~isempty(yinit))             % An initial configuration exists.
        y = yinit;
      elseif (strcmp(opt.init, 'cs')) & (~isempty(winit))
        y = D*winit;
      else   
        y = mds_init(D,n,opt.init);
      end

      if (~isempty(opt.nanid))         % Rescale.
        scale = (max(yrep)-min(yrep))./(max(y)-min(y));
        y = y .* repmat(scale,ni,1); 
      end

      % The final number of true distances is:
      no_dist = (ni*nc - length(opt.nanid));
    end

    % Check whether there is enough data left.

    if (~isempty(opt.nanid))
      if (n*ni+2 > no_dist),
        error('There are too many missing distances: it is not possible to determine the MDS map.');
      end
      if (strcmp (opt.itmap,'no'))
        opt.itmap = 'yes';
        prwarning(1,'Due to the missing values, the projection can only be iterative. OPTIONS are changed appropriately.')
      end
    end

    if (opt.inspect > 0)
      opt.plotlab = lab(I,:);  % Labels to be used for plotting in MDS_SAMMAP.
    else
      opt.plotlab = [];
    end                       

    if (isempty(yrep)) | (~isempty(yrep) & strcmp(opt.itmap, 'yes'))  

      % Either no MDS exists yet OR new points should be mapped in an iterative manner.

      printinfo(opt);
      [yy,err] = mds_sammap(D,y,yrep,opt);

      % Define the linear projection of distances.

      v = [];
      if (isempty(yrep)) & (isempty(opt.nanid))  
        if (rank(D) < m)
          v = prpinv(D)*yy;
        else
          v = D \ yy;
        end
      end
    else
      % New points should be added by a linear projection of dissimilarity data.
      yy = D*v;
    end

    % Establish the projected configuration including the removed points.

    y = zeros(m,n); y(I,:) = +yy;                   
    if (~isempty(J))
      if (~isempty(yrep))
        y(J,:) = +yrep(II(P),:);
      else
        for k=length(J):-1:1        % J: indices of removed points.
          y(J(k),:) = y(P(k),:);    % P: indices of points equal to points in J.
        end 
      end
    end

    % In the definition step: shift the obtained configuration such that the
    % mean lies at the origin.

    if (isempty(yrep))
      y = y - ones(length(y),1)*mean(y);
      y = prdataset(y,lab);
    else
      w = prdataset(y,lab);
      return;
    end

    % In the definition step: the mapping should be stored.

    opt = rmfield(opt,'nanid');   % These fields are to be used internally only;
    opt = rmfield(opt,'plotlab'); % not to be set up from outside
    w   = prmapping(mfilename,'trained',{y,opt,I,winit,v,[]},[],m,n);
    w   = setname(w,mapname);
    
  end

return

% **********************************************************************************
%                             Extra functions
% **********************************************************************************

% PRINTINFO(OPT)
%
% Prints progress information to the file handle given by OPT.ST.

function printinfo(opt)

  if opt.st < 1
    return
  end
  %fprintf(opt.st,'Sammon mapping, error function with the parameter q=%d\n',opt.q);

  switch (opt.optim)
    case 'pn',  
      %fprintf(opt.st,'Minimization by Pseudo-Newton algorithm\n');
    case 'scg',  
      %fprintf(opt.st,'Minimization by Scaled Conjugate Gradients algorithm\n');
    otherwise 
      error(strcat('Possible initialization methods: pn (Pseudo-Newton), ',...
                   'or scg (Scaled Conjugate Gradients).'));
  end

return

% **********************************************************************************

% [I,J,P] = MDS_REPPOINTS(D)
%
% Finds the indices of repeated/left points. J contains the indices of
% repeated points in D. This means that for each j in J, there is a point
% k ~= j such that D(J(j),k) = 0. I contains the indices of the remaining 
% points, and P those of the points in I that were identical to those in J.
% Directly used in MDS routine.

function [I,J,P] = mds_reppoints(D)

  epsilon = 1e-20;      % Differences smaller than this are assumed to be zero.

  [m,mm] = size(D);

  I = 1:m; J = []; P = [];

  if (m == mm) & (all(abs(D(1:m+1:end)) <= epsilon))
    K = intersect (find (triu(ones(m),1)), find(D < epsilon));
    if (~isempty(K))
      P = mod(K,m);
      J = fix(K./m) + 1;           
      I(J) = [];                    
    end
  else
    [J,P] = find(D<=epsilon); 
    I(J)  = [];
  end

return;

% **********************************************************************************

% MDS_SETOPT  Sets parameters for the MDS mapping
%
%   OPT = MDS_SETOPT(OPT_GIVEN)
% 
% INPUT
%   OPT_GIVEN Parameters for the MDS mapping, described below (default: [])
%
% OUTPUT
%   OPT       Structure of chosen options; if OPT_GIVEN is empty, then 
%               OPT.q       = 0
%               OPT.optim   = 'pn'
%               OPT.init    = 'cs'
%               OPT.etol    = 1e-6 (the precise value depends on q)
%               OPT.maxiter = 100 
%               OPT.itmap   = 'yes' 
%               OPT.isratio = 0 
%               OPT.st      = 1 
%               OPT.inspect = 2 
%
% DESCRIPTION  
% Parameters for the MDS mapping can be set or changed. OPT_GIVEN consists
% of the following fields: 'q', 'init', 'optim', 'etol', 'maxiter','itmap',
% 'isratio', 'st' and 'inspect'. OPTIONS can include all the fields or some
% of them only. The fields of OPT have some default values, which can be
% changed by the OPT_GIVEN field values. If OPT_GIVEN is empty, then OPT
% contains all default values. For a description of the fields, see MDS.

%
% Copyright: Elzbieta Pekalska, ela@ph.tn.tudelft.nl, 2000-2003
% Faculty of Applied Sciences, Delft University of Technology
%

function opt = mds_setopt (opt_given)

  opt.q       = 0;
  opt.init    = 'cs'; 
  opt.optim   = 'pn';
  opt.st      = 1; 
  opt.itmap   = 'yes'; 
  opt.maxiter = 100; 
  opt.inspect = 2; 
  opt.etol    = inf; 
  opt.isratio = 0; 
  opt.nanid   = [];    % Here are some extra values; set up in the MDS routine.
  opt.plotlab = [];    % Not to be changed by the user.

  if (~isempty(opt_given))
    if (~isstruct(opt_given))
      error('OPTIONS should be a structure with at least one of the following fields: q, init, etol, optim, maxiter, itmap, isratio, st or inspect.');
    end
    fn = fieldnames(opt_given);
    if (~all(ismember(fn,fieldnames(opt))))
      error('Wrong field names; valid field names are: q, init, optim, etol, maxiter, itmap, isratio, st or inspect.')
    end 
    for i = 1:length(fn)
      opt = setfield(opt,fn{i},getfield(opt_given,fn{i}));
    end
  end

  if (isempty(intersect(opt.q,-2:1:2)))
    error ('OPTIONS.q should be -2, -1, 0, 1 or 2.');
  end

  if (opt.maxiter < 2)
    error ('OPTIONS.iter should be at least 1.');
  end

  if (isinf(opt.etol))
    switch (opt.q)                 % Different defaults for different stresses.
      case -2, 
        opt.etol = 1e-6;
      case {-1,0}
        opt.etol = 10*sqrt(eps); 
      case {1,2}
        opt.etol = sqrt(eps); 
      end 
  elseif (opt.etol <= 0) | (opt.etol >= 0.1)
    error ('OPTIONS.etol should be positive and smaller than 0.1.');
  end

  if (~ismember(opt.optim, {'pn','scg'}))
    error('OPTIONS.optim should be pn or scg.');
  end

  if (~ismember(opt.itmap, {'yes','no'}))
    error('OPTIONS.itmap should be yes or no.');
  end

return


% MDS_SAMMAP Sammon iterative nonlinear mapping for MDS
%
%   [YMAP,ERR] = MDS_SAMMAP(D,Y,YREP,OPTIONS)
% 
% INPUT
%   D       Square (M x M) dissimilarity matrix
%   Y       M x N matrix containing starting configuration, or
%   YREP    Configuration of the representation points
%   OPTIONS Various parameters of the  minimization procedure put into 
%           a structure consisting of the following fields: 'q', 'optim',
%           'init','etol','maxiter', 'isratio', 'itmap', 'st' and 'inspect'
%           (default:
%             OPTIONS.q       = 0          
%             OPTIONS.optim   = 'pn'      
%             OPTIONS.init    = 'cs'
%             OPTIONS.etol    = 1e-6 (the precise value depends on q)
%             OPTIONS.maxiter = 100 
%             OPTIONS.isratio = 0 
%             OPTIONS.itmap   = 'yes' 
%             OPTIONS.st      = 1 
%             OPTIONS.inspect = 2).
% 
% OUTPUT
%   YMAP     Mapped configuration
%   ERR     Sammon stress 
%
% DESCRIPTION  
% Maps the objects given by a symmetric distance matrix D (with a zero
% diagonal) onto, say, an N-dimensional configuration YMAP by an iterative
% minimization of a variant of the Sammon stress. The minimization starts
% from the initial configuration Y; see MDS_INIT.
%
% YREP is the Sammon configuration of the representation set. It is used
% when new points have to be projected. In other words, if D is an M x M
% symmetric distance matrix, then YREP is empty; if D is an M x N matrix,
% then YMAP is sought such that D can approximate the distances between YMAP
% and YREP.
%
% Missing values can be handled by marking them by NaN in D.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, MDS, MDS_CS, MDS_INIT, MDS_SETOPT

%
% Copyright: Elzbieta Pekalska, ela@ph.tn.tudelft.nl, 2000-2003
% Faculty of Applied Sciences, Delft University of Technology
%

function [y,err] = mds_sammap(Ds,y,yrep,opt)

  if (nargin < 4)
    opt = [];                 % Will be filled by MDS_SETOPT, below.
  end
  if isempty(opt), opt = mds_setopt(opt); end
  if (nargin < 3)
    yrep = []; 
  end

  % Extract labels and calculate distance matrix.

  [m,n] = size(y);
  if (~isempty(yrep))
    replab = getlab(yrep);
    yrep = +yrep;
    D = sqrt(distm(y,yrep)); 
  else
    D = sqrt(distm(y));     
    yrep = [];
    replab = [];
  end

%   if (isempty(opt.plotlab))
%     opt.plotlab = ones(m,1);
%   end

  it   = 0;           % Iteration number.
  eold = inf;         % Previous stress (error).

  % Calculate initial stress.

  [e,a]= mds_samstress(opt.q,y,yrep,Ds,D,opt.nanid,opt.isratio); err = e;

  %if opt.st > 0, fprintf(opt.st,'iteration: %4i   stress: %3.8d\n',it,e); end
  
  tt = sprintf('Processing %i iterations: ',opt.maxiter);
  prwaitbar(opt.maxiter,tt)

  switch (opt.optim)  
    case 'pn',        % Pseudo-Newton minimization.

      epr   = e;      % Previous values of E, EOLD for line search algorithm.
      eepr  = e;
      add   = 1e-4;   % Avoid G2 equal 0; see below.
      
      BETA   = 1e-3;   % Parameter for line search algorithm.
      lam1   = 0;     
      lam2   = 1;      % LAM (the step factor) lies in (LAM1, LAM2).
      LMIN   = 1e-10;  % Minimum acceptable value for LAM in line search.
      lambest = 1e-4; % LAM = LAMBEST, if LAM was not found.

      % Loop until the error change falls below the tolerance or until
      % the maximum number of iterations is reached.
      
      while (abs(e-eold) >= opt.etol*(1 + abs(e))) & (it < opt.maxiter)

        % Plot progress if requested.

%        if (opt.st > 0) & (opt.inspect > 0) & (mod(it,opt.inspect) == 0)
%          if (it == 0), figure(1); clf; end
%          mds_plot(y,yrep,opt.plotlab,replab,e); 
%        end

        yold = y; eold = e;

        % Calculate the stress.

        [e,a]      = mds_samstress (opt.q,yold,yrep,Ds,[],opt.nanid,opt.isratio);

        % Calculate gradients and pseudo-Newton update P.

        [g1,g2,cc] = mds_gradients (opt.q,yold,yrep,a*Ds,[],opt.nanid);
        p        = (g1/cc)./(abs(g2/cc) + add);
        slope    = g1(:)' * p(:);

        lam = 1;                 

        % Search for a suitable step LAM using a line search.

        while (1)

          % Take a step and calculate the delta error DE.

          y      = yold + lam .* p;  
          [e,a]  = mds_samstress(opt.q,y,yrep,Ds,[],opt.nanid,opt.isratio);
          de     = e - eold;

          % Stop if the condition for a suitable lam is fulfilled.

          if (de < BETA * lam * slope)
            break;
          end
   
          % Try to find a suitable step LAM.

          if (lam ~= 1)
            r1 = de - lam*slope; 
            r2 = epr - eepr - lam2*slope; 
            aa = (2*de + r1/lam^2 - r2/lam2^2)/(lam2-lam1)^3;
            bb = ((-lam2*r1)/lam^2 +(lam*r2)/lam2^2)/(lam2-lam1)^2;

            if (abs(aa) <= eps)
              lamtmp = -0.5 * slope/bb;
            else
              lamtmp = (-bb + sqrt(max(0,(bb^2 - 3*aa*slope))))/(3*aa);  
            end

            % Prevent LAM from becoming too large.

            if (lamtmp > 0.5 * lam), lamtmp = 0.5 * lam; end

          else
            lamtmp = -0.5 * slope/(e - eold - slope);
          end

          % Prevent LAM from becoming too small.

          lam2 = lam;    
          lam  = max(lamtmp, 0.1*lam);  
          epr  = e;
          eepr = eold;   

          if (lam < LMIN)
            y = yold + lambest .* p;
            [e,a] = mds_samstress(opt.q,y,yrep,Ds,[],opt.nanid,opt.isratio);
            break;
          end
        end

        it  = it + 1; 
        err = [err; e];
        
        prwaitbar(opt.maxiter,it,[tt num2str(it) ', error: ' num2str(e)]);
        %if opt.st > 0, fprintf(opt.st,'iteration: %4i   stress: %3.8d\n',it,e); end

      end

    case 'scg',              % Scaled Conjugate Gradient minimization.

      sigma0  = 1e-8;                        
      lambda  = 1e-8;       % Regularization parameter.
      lambda1 = 0;                                                             

      % Calculate initial stress and direction of decrease P.

      [e,a]   = mds_samstress(opt.q,y,yrep,Ds,D,opt.nanid,opt.isratio);   
      g1      = mds_gradients(opt.q,y,yrep,a*Ds,D,opt.nanid);   % Gradient
      p       = -g1;                                            % Direction of decrease
      gnorm2  = g1(:)' * g1(:);
      success = 1;                                              

      % Sanity check.

      if (gnorm2 < 1e-15)
        prwarning(2,['Gradient is nearly zero: ' gnorm2 ' (unlikely); initial configuration is the right one.']);  
        return
      end

      % Loop until the error change falls below the tolerance or until
      % the maximum number of iterations is reached.

      while (abs(eold - e) > opt.etol * (1 + e)) & (it < opt.maxiter)

        g0     = g1;            % Previous gradient
        pnorm2 = p(:)' * p(:);
        eold   = e;

        % Plot progress if requested.

%        if (opt.inspect > 0) & (mod(it,opt.inspect) == 0) 
%          if (it == 0), figure(1); clf; end
%          mds_plot(y,yrep,opt.plotlab,replab,e); 
%        end

        if (success)
          sigma = sigma0/sqrt(pnorm2);      % SIGMA: a small step from y to yy
          yy    = y + sigma .*p;
          [e,a] = mds_samstress(opt.q,yy,yrep,Ds,[],opt.nanid,opt.isratio); 
          g2    = mds_gradients(opt.q,yy,yrep,a*Ds,[],opt.nanid); % Gradient for yy
          s     = (g2-g1)/sigma;        % Approximation of Hessian*P directly, instead of computing
          delta = p(:)' * s(:);         % the Hessian, since only DELTA = P'*Hessian*P is in fact needed
        end

        % Regularize the Hessian indirectly to make it positive definite.
        % DELTA is now computed as P'* (Hessian + regularization * Identity) * P
        delta = delta + (lambda1 - lambda) * pnorm2;  

        % Indicate if the Hessian is negative definite; the regularization above was too small.
        if (delta < 0)
          lambda1 = 2 * (lambda - delta/pnorm2);  % Now the Hessian will be positive definite
          delta   = -delta + lambda * pnorm2;     % This is obtained after plugging lambda1 
          lambda  = lambda1;                      % into the formulation a few lines above.
        end

        mi = - p(:)' * g1(:);
        yy = y + (mi/delta) .*p;                  % mi/delta is a step size 
        [ee,a] = mds_samstress(opt.q,yy,yrep,Ds,[],opt.nanid,opt.isratio); 

        % This minimization procedure is based on the second order approximation of the  
        % stress function by using the gradient and the Hessian approximation. The Hessian 
        % is regularized, but maybe not sufficiently. The ratio Dc (<=1) below indicates 
        % a proper approximation, if Dc is close to 1.

        Dc = 2 * delta/mi^2 * (e - ee); 
        e  = ee;

        % If Dc > 0, then the stress can be successfully decreased.
        success = (Dc >= 0);
        if (success)
          y       = yy;
          [ee,a]  = mds_samstress(opt.q,yy,yrep,Ds,[],opt.nanid,opt.isratio); 
          g1      = mds_gradients(opt.q,yy,yrep,a*Ds,[],opt.nanid); 
          gnorm2  = g1(:)' * g1(:);
          lambda1 = 0;

          beta = max(0,(g1(:)'*(g1(:)-g0(:)))/mi);  
          p   = -g1 + beta .* p;                 % P is a new conjugate direction  

          if (g1(:)'*p(:) >= 0 | mod(it-1,n*m) == 0),
            p = -g1;                            % No much improvement, restart
          end   

          if (Dc >= 0.75)
            lambda = 0.25 * lambda;             % Make the regularization smaller
          end

          it  = it + 1; 
          err = [err; e];
          prwaitbar(opt.maxiter,it,[tt num2str(it) ', error: ' num2str(e)]);
          %fprintf (opt.st,'iteration: %4i   stress: %3.8d\n',it,e); 

        else        % Dc < 0  
          % Note that for Dc < 0, the iteration number IT is not increased, 
          % so the Hessian is further regularized until SUCCESS equals 1.
          lambda1 = lambda;
        end

      % The approximation of the Hessian was poor or the stress was not 
      % decreased (Dc < 0), hence the regularization lambda is enlarged.
      if (Dc < 0.25)
        lambda = lambda + delta * (1 - Dc)/pnorm2;
      end
    end
  end
  prwaitbar(0);

return


% **********************************************************************************

%MDS_STRESS - Calculate Sammon stress during optimization
%
%   E = MDS_SAMSTRESS(Q,Y,YREP,Ds,D,NANINDEX)
%
% INPUT
%   Q         Indicator of the Sammon stress; Q = -2,-1,0,1,2
%   Y         Current lower-dimensional configuration
%   YREP      Configuration of the representation objects; it should be
%             empty when no representation set is considered
%   Ds        Original distance matrix
%   D         Approximate distance matrix (optional; otherwise computed from Y)
%   NANINDEX  Index of the missing values; marked in Ds by NaN (optional; to
%             be found in Ds)
%
% OUTPUT
%   E         Sammon stress
%
% DESCRIPTION
% Computes the Sammon stress between the original distance matrix Ds and the
% approximated distance matrix D between the mapped configuration Y and the
% configuration of the representation set YREP, expressed as follows:
%
% E = 1/(sum_{i<j} Ds_{ij}^(q+2)) sum_{i<j} (Ds_{ij} - D_{ij})^2 * Ds_{ij}^q
%
% It is directly used in the MDS_SAMMAP routine.

%
% Copyright: Elzbieta Pekalska, Robert P.W. Duin, ela@ph.tn.tudelft.nl, 2000-2003
% Faculty of Applied Sciences, Delft University of Technology
%

function [e,alpha] = mds_samstress (q,y,yrep,Ds,D,nanindex,isratio)

  % If D is not given, calculate it from Y, the current mapped points.

  if (nargin < 5) | (isempty(D))
    if (~isempty(yrep))
      D = sqrt(distm(y,yrep));     
    else
      D = sqrt(distm(y));     
    end
  end
  
  if (nargin < 6)
    nanindex = [];    % Not given, so calculate below.
  end

  if (nargin < 7)
    isratio = 0;      % Assume this is meant.
  end

  todefine = isempty(yrep);

  [m,n]  = size(y); [mm,k] = size(Ds);
  if (m ~= mm)
    error('The sizes of Y and Ds do not match.');
  end
  if (any(size(D) ~= size(Ds)))
    error ('The sizes of D and Ds do not match.');
  end
  m2 = m*k;

  % Convert to double.
  D  = +D; Ds = +Ds;

  % I is the index of non-NaN, non-zero (> eps) values to be included 
  % for the computation of the stress.

  I = 1:m2; 
  if (~isempty(nanindex))
    I(nanindex) = [];
  end

  O = [];
  if (todefine), O = 1:m+1:m2; end                               
  I = setdiff(I,O);

  % If OPTIONS.isratio is set, calculate optimal ALPHA to scale with.

  if (isratio)
    alpha = sum((Ds(I).^q).*D(I).^2)/sum((Ds(I).^(q+1)).*D(I));
    Ds    = alpha*Ds;
  else
    alpha = 1; 
  end
    
  % C is the normalization factor.
  c = sum(Ds(I).^(q+2)); 

  % If Q == 0, prevent unnecessary calculation (X^0 == 1).
  if (q ~= 0)
    e = sum(Ds(I).^q .*((Ds(I)-D(I)).^2))/c;
  else
    e = sum(((Ds(I)-D(I)).^2))/c;
  end

return


% **********************************************************************************

%MDS_GRADIENTS - Gradients for variants of the Sammon stress
%
%   [G1,G2,CC] = MDS_GRADIENTS(Q,Y,YREP,Ds,D,NANINDEX)
%
% INPUT
%   Q          Indicator of the Sammon stress; Q = -2,-1,0,1,2
%   Y         Current lower-dimensional configuration
%   YREP      Configuration of the representation objects; it should be
%             empty when no representation set is considered
%   Ds        Original distance matrix
%   D         Approximate distance matrix (optional; otherwise computed from Y)
%   nanindex  Index of missing values; marked in Ds by NaN (optional;
%             otherwise found from Ds)
%
% OUTPUT
%   G1        Gradient direction
%   G2        Approximation of the Hessian by its diagonal  
%
% DESCRIPTION  
% This is a routine used directly in the MDS_SAMMAP routine.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MDS, MDS_INIT, MDS_SAMMAP, MDS_SAMSTRESS

%
% Copyright: Elzbieta Pekalska, Robert P.W. Duin, ela@ph.tn.tudelft.nl, 2000-2003
% Faculty of Applied Sciences, Delft University of Technology
%

function [g1,g2,c] = mds_gradients(q,y,yrep,Ds,D,nanindex)

  % If D is not given, calculate it from Y, the current mapped points.

	y = +y;
	yrep = +yrep;
  if (nargin < 5) | (isempty(D))
    if (~isempty(yrep))
      D = sqrt(distm(y,yrep));     
    else
      D = sqrt(distm(y));     
    end
  end

  % If NANINDEX is not given, find it from Ds.

  if (nargin < 6)
    nanindex = find(isnan(Ds(:))==1);  
  end

  % YREP is empty if no representation set is defined yet.
  % This happens when the mapping should be defined.

  todefine = (isempty(yrep));
  if (todefine) 
    yrep  = y;
  end

  [m,n]  = size(y); [mm,k] = size(Ds);
  if (m ~= mm)
    error('The sizes of Y and Ds do not match.');
  end
  if (any(size(D) ~= size(Ds)))
    error ('The sizes of D and Ds do not match.');
  end
  m2 = m*k;

  % Convert to doubles.

  D  = +D; Ds = +Ds;

  % I is the index of non-NaN, non-zero (> eps) values to be included 
  % for the computation of the gradient and the Hessians diagonal.

  I = (1:m2)';
  if (~isempty(nanindex))
    I(nanindex) = [];
  end
  K    = find(Ds(I) <= eps | D(I) <= eps);
  I(K) = [];

  % C is a normalization factor.

  c = -2/sum(Ds(I).^(q+2)); 

  % Compute G1, the gradient.

  h1 = zeros(m2,1);
  if (q == 0)                    % Prevent unnecessary computation when Q = 0.
    h1(I) = (Ds(I)-D(I)) ./ D(I);  
  else  
    h1(I) = (Ds(I)-D(I)) ./ (D(I).*(Ds(I).^(-q)));
  end
  h1 = reshape (h1',m,k);
  g2 = h1 * ones(k,n);          % Here G2 is assigned only temporarily,
  g1 = c * (g2.*y - h1*yrep);    % for the computation of G1.     

  % Compute G2, the diagonal of the Hessian, if requested.

  if (nargout > 1)
    h2 = zeros(m2,1);
    switch (q)
      case -2, 
        h2(I) = -1./(Ds(I).*D(I).^3);
      case -1, 
        h2(I) = -1./(D(I).^3);
      case 0, 
        h2(I) = - Ds(I)./(D(I).^3);
      case 1, 
        h2(I) = - Ds(I).^2./(D(I).^3);
      case 2, 
        h2(I) = -(Ds(I)./D(I)).^3;
      end 
    h2 = reshape (h2',m,k);
    g2 = c * (g2 + (h2*ones(k,n)).*y.^2 + h2*yrep.^2 - 2*(h2*yrep).*y);  
  end

return
