%OUT2 Mapping that returns second output of input mapping
%
%   V = OUT2(W)
%   V = W*OUT2
%
% INPUT
%   W   Fixed or trained mapping
%
% OUTPUT
%   V   Fixed out trained mapping
%
% DESCRIPTION
% In case A*W (A is a dataset or double) returns two outputs then
% A*V returns just the second output of these two.
%
% EXAMPLE
% test2 = testc*out2;% define testc for second output par (# class errors)
% a = gendatd;       % train set
% t = gendatd;       % test set
% t*knnc(a,1)*test2  % execute, list of # errors per class
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function out = out2(par1,par2)

if nargin == 0
  out = define_mapping([],'combiner');
elseif nargin == 1 & ismapping(par1)
  out = define_mapping({[],par1},'fixed',getname(par1));
elseif nargin == 2
  [dummy,out] = par1*par2;
else
  error('Illegal input');
end