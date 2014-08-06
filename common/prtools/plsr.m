% PLSR Partial Least Squares Regression
%
%  W = PLSR
%  W = PLSR([],MAXLV,METHOD)
%
%  [W, INFORM] = PLSR(A,MAXLV,METHOD)
%
% INPUT 
%   A             training dataset
%   MAXLV         maximal number of latent variables (will be corrected
%                 if > rank(A)); 
%                 MAXLV=inf means MAXLV=min(size(A)) -- theoretical
%                 maximum number of LV; 
%                 by default = inf
%   METHOD       'NIPALS' or 'SIMPLS'; by default = 'SIMPLS'
%
% OUTPUT 
%   W             PLS feature extraction mapping 
%   INFORM        extra algorithm output
%
% DESRIPTION
% PRTools Adaptation of PLS_TRAIN/PLS_APPLY routines. No preprocessing
% is done inside this mapping. It is the user responsibility to train
% preprocessing on training data and apply it to the test data.
%
% Crisp labels will be converted into soft labels which will be used as
% a target matrix.
%
% In order to do regression with the smaller number of latent variables
% than the number of LV's selected during trainig do
% d = w.data;
% d.n = new_n;
% w.data = d;
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PLS_TRAIN, PLS_TRANSFORM, PLS_APPLY

% Copyright: S.Verzakov, s.verzakov@ewi.tudelft.nl 
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: plsr.m,v 1.1 2007/08/28 11:00:39 davidt Exp $
%

function [w,inform]=plsm(par1,par2,par3)
	
% No dataset given: return untrained mapping.
if (nargin < 1) | (isempty(par1))
  if nargin < 2
    par2 = inf;
  end  
  if nargin < 3
    par3 = 'SIMPLS';
  end  
  data = {par2,par3};
  w = prmapping(mfilename,'untrained',data);
	w = setname(w,'Partial Least Squares Regression');
	return
end

%isdataset(par1);           % Assert that A is a dataset.

% training
if nargin < 2 | ~isa(par2,'prmapping') 
  % a*w when w is untrained or 
  if nargin < 2
    par2 = inf;
  end  
  if nargin < 3
    par3 = 'SIMPLS';
  end  
  maxLV  = par2;
  method = par3;
 
	if strcmp(par1.labtype,'crisp')
	  y=gettargets(setlabtype(par1,'soft'));
	else
  	y=gettargets(par1);
  end

	% options
  Options.maxLV  = maxLV;
  Options.method = method;
	Options.X_centering=[];
  Options.Y_centering=[];
	Options.X_scaling=[];
  Options.Y_scaling=[];

	[B,XRes,YRes,Options]=pls_train(+par1,y,Options);
	
  clear B

  data.n=Options.maxLV;
  data.R=XRes.R;
  data.C=YRes.C;
  data.Options=Options;
	
	% Save all useful data.
	w = prmapping(mfilename,'trained',data,[],size(XRes.R,1),size(YRes.C,1));
	w = setname(w,'Partial Least Squares Mapping');

	if nargout > 1
	 inform.XRes=XRes;
	 inform.YRes=YRes;
	end
	
% execution
else 
	data = getdata(par2); % Unpack the mapping.
    if data.n > size(data.R,2)
    ErrMsg = sprintf(['PLS: The nubmer of LV(s) asked (%d) is greater than\n'...
                      'the number of LV(s) available (%d)'],data.n,size(data.R,2));
    error(ErrMsg);
  end

  Y = pls_apply(+par1,data.R(:,1:data.n)*data.C(:,1:data.n)',data.Options);
  w = setdat(par1,Y,par2);
end

return

