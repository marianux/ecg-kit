function out1 = pr_newrb(varargin)
%
%  A copy of Matlab's newrb, needed to correct a bug:
%  displayFreq = inf is not accepted, by which display could not be skipped
%
%NEWRB Design a radial basis network.
%
%  Radial basis networks can be used to approximate functions.  <a href="matlab:doc newrb">newrb</a>
%  adds neurons to the hidden layer of a radial basis network until it
%  meets the specified mean squared error goal.
%
%  <a href="matlab:doc newrb">newrb</a>(X,T,GOAL,SPREAD,MN,DF) takes these arguments,
%    X      - RxQ matrix of Q input vectors.
%    T      - SxQ matrix of Q target class vectors.
%    GOAL   - Mean squared error goal, default = 0.0.
%    SPREAD - Spread of radial basis functions, default = 1.0.
%    MN     - Maximum number of neurons, default is Q.
%    DF     - Number of neurons to add between displays, default = 25.
%  and returns a new radial basis network.
%
%  The larger that SPREAD is the smoother the function approximation
%  will be.  Too large a spread means a lot of neurons will be
%  required to fit a fast changing function.  Too small a spread
%  means many neurons will be required to fit a smooth function,
%  and the network may not generalize well.  Call NEWRB with
%  different spreads to find the best value for a given problem.
%
%  Here we design a radial basis network given inputs X and targets T.
%
%    X = [1 2 3];
%    T = [2.0 4.1 5.9];
%    net = <a href="matlab:doc newrb">newrb</a>(X,T);
%    Y = net(X)
%
%  See also SIM, NEWRBE, NEWGRNN, NEWPNN.

% Mark Beale, 11-31-97
% Copyright 1992-2010 The MathWorks, Inc.
% $Revision: 1.1.6.9 $ $Date: 2010/04/24 18:10:03 $

%% =======================================================
%  BOILERPLATE_START
%  This code is the same for all Network Functions.

  persistent INFO;
  if isempty(INFO), INFO = get_info; end
  if (nargin > 0) && ischar(varargin{1})
    code = varargin{1};
    switch code
      case 'info',
        out1 = INFO;
      case 'check_param'
        err = check_param(varargin{2});
        if ~isempty(err), nnerr.throw('Args',err); end
        out1 = err;
      case 'create'
        if nargin < 2, nnerr.throw('Not enough arguments.'); end
        param = varargin{2};
        err = nntest.param(INFO.parameters,param);
        if ~isempty(err), nnerr.throw('Args',err); end
        out1 = create_network(param);
        out1.name = INFO.name;
      otherwise,
        % Quick info field access
        try
          out1 = eval(['INFO.' code]);
        catch %#ok<CTCH>
          nnerr.throw(['Unrecognized argument: ''' code ''''])
        end
    end
  else
    [param,err] = INFO.parameterStructure(varargin);
    if ~isempty(err), nnerr.throw('Args',err); end
    net = create_network(param);
    net.name = INFO.name;
    out1 = init(net);
  end
end

function v = fcnversion
  v = 7;
end

%  BOILERPLATE_END
%% =======================================================

function info = get_info
  info = nnfcnNetwork(mfilename,'Radial Basis Network',fcnversion, ...
    [ ...
    nnetParamInfo('inputs','Input Data','nntype.data',{0},...
    'Input data.'), ...
    nnetParamInfo('targets','Target Data','nntype.data',{0},...
    'Target output data.'), ...
    nnetParamInfo('goal','Performance Goal','nntype.pos_scalar',0,...
    'Performance goal.'), ...
    nnetParamInfo('spread','Radial basis spread','nntype.strict_pos_scalar',1,...
    'Distance from radial basis center to 0.5 output.'), ...
    nnetParamInfo('maxNeurons','Maximum number of neurons','nntype.pos_int_inf_scalar',inf,... % TODO - type
    'Maximum number of neurons to add to network.'), ...
    nnetParamInfo('displayFreq','Display Frequency','nntype.pos_int_inf_scalar',1,...
    'Number of added neurons between displaying progress at command line.'), ...
    ]);
end

function err = check_param(param)
  err = '';
end

function net = create_network(param)

  % Data
  p = param.inputs;
  t = param.targets;
  if iscell(p), p = cell2mat(p); end
  if iscell(t), t = cell2mat(t); end

  % Max Neurons
  Q = size(p,2);
  mn = param.maxNeurons;
  if (mn > Q), mn = Q; end

  % Dimensions
  R = size(p,1);
  S2 = size(t,1);

  % Architecture
  net = network(1,2,[1;1],[1; 0],[0 0;1 0],[0 1]);

  % Simulation
  net.inputs{1}.size = R;
  net.layers{1}.size = 0;
  net.inputWeights{1,1}.weightFcn = 'dist';
  net.layers{1}.netInputFcn = 'netprod';
  net.layers{1}.transferFcn = 'radbas';
  net.layers{2}.size = S2;
  net.outputs{2}.exampleOutput = t;

  % Performance
  net.performFcn = 'mse';

  % Design Weights and Bias Values
  warn1 = warning('off','MATLAB:rankDeficientMatrix');
  warn2 = warning('off','MATLAB:nearlySingularMatrix');
  [w1,b1,w2,b2,tr] = designrb(p,t,param.goal,param.spread,mn,param.displayFreq);
  warning(warn1.state,warn1.identifier);
  warning(warn2.state,warn2.identifier);

  net.layers{1}.size = length(b1);
  net.b{1} = b1;
  net.iw{1,1} = w1;
  net.b{2} = b2;
  net.lw{2,1} = w2;
end

%======================================================
function [w1,b1,w2,b2,tr] = designrb(p,t,eg,sp,mn,df)

  [r,q] = size(p);
  [s2,q] = size(t);
  b = sqrt(-log(.5))/sp;

  % RADIAL BASIS LAYER OUTPUTS
  P = radbas(dist(p',p)*b);
  PP = sum(P.*P)';
  d = t';
  dd = sum(d.*d)';

  % CALCULATE "ERRORS" ASSOCIATED WITH VECTORS
  e = ((P' * d)' .^ 2) ./ (dd * PP');

  % PICK VECTOR WITH MOST "ERROR"
  pick = findLargeColumn(e);
  used = [];
  left = 1:q;
  W = P(:,pick);
  P(:,pick) = []; PP(pick,:) = [];
  e(:,pick) = [];
  used = [used left(pick)];
  left(pick) = [];

  % CALCULATE ACTUAL ERROR
  w1 = p(:,used)';
  a1 = radbas(dist(w1,p)*b);
  [w2,b2] = solvelin2(a1,t);
  a2 = w2*a1 + b2*ones(1,q);
  MSE = mse(t-a2);

  % Start
  tr = nntraining.newtr(mn,'perf');
  tr.perf(1) = mse(t-repmat(mean(t,2),1,q));
  tr.perf(2) = MSE;
  if isfinite(df)
    fprintf('NEWRB, neurons = 0, MSE = %g\n',tr.perf(1));
  end
  flag_stop = 0;

  iterations = min(mn,q);
  for k = 2:iterations

    % CALCULATE "ERRORS" ASSOCIATED WITH VECTORS
    wj = W(:,k-1);
    a = wj' * P / (wj'*wj);
    P = P - wj * a;
    PP = sum(P.*P)';
    e = ((P' * d)' .^ 2) ./ (dd * PP');

    % PICK VECTOR WITH MOST "ERROR"
    pick = findLargeColumn(e);
    W = [W, P(:,pick)];
    P(:,pick) = []; PP(pick,:) = [];
    e(:,pick) = [];
    used = [used left(pick)];
    left(pick) = [];

    % CALCULATE ACTUAL ERROR
    w1 = p(:,used)';
    a1 = radbas(dist(w1,p)*b);
    [w2,b2] = solvelin2(a1,t);
    a2 = w2*a1 + b2*ones(1,q);
    MSE = mse(t-a2);

    % PROGRESS
    tr.perf(k+1) = MSE;

    % DISPLAY
    if isfinite(df) & (~rem(k,df))
      fprintf('NEWRB, neurons = %g, MSE = %g\n',k,MSE);
      flag_stop=plotperf(tr,eg,'NEWRB',k);
    end

    % CHECK ERROR
    if (MSE < eg), break, end
    if (flag_stop), break, end

  end

  [S1,R] = size(w1);
  b1 = ones(S1,1)*b;

  % Finish
  tr = nntraining.cliptr(tr,k);
end

%======================================================
function i = findLargeColumn(m)
  replace = find(isnan(m));
  m(replace) = zeros(size(replace));
  m = sum(m .^ 2,1);
  i = find(m == max(m));
  i = i(1);
end

%======================================================

function [w,b] = solvelin2(p,t)
  if nargout <= 1
    w= t/p;
  else
    [pr,pc] = size(p);
    x = t/[p; ones(1,pc)];
    w = x(:,1:pr);
    b = x(:,pr+1);
  end
end

%======================================================
