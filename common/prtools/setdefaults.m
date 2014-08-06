%SETDEFAULTS Set defaults for function parameters
%
%   ARGOUT = SETDEFAULTS(ARGIN,DEF1,DEF2, ....)
%   [P1,P2, ...] = SETDEFAULTS(ARGIN,DEF1,DEF2, ....)
%
% INPUT
%   ARGIN   - Cell array with function input arguments, typically VARARGIN
%   DEF1    - Default value for argument 1
%   DEF2    - Default value for argument 2
%
% OUTPUT
%   ARGOUT  - Cell array defaults replacing the empty input arguments
%   P1      - Input argument 1, replaced by its default if empty
%   P2      - Input argument 2, replaced by its default if empty
%
% DESCRIPTION
% This routine substitutes empty input parameters of a function (typically
% given by VARARGIN) with the defaults DEF1, DEF2, etcetera.

% Copyright: Robert P.W. Duin, prtools@rduin.nl

function varargout = setdefaults(parin,varargin)

if ~iscell(parin), parin = {parin}; end
varin = cell(1,max([nargout,numel(parin),numel(varargin)]));
varin(1:numel(varargin)) = varargin;

varout = cell(1,numel(varin));
varout(1:numel(parin)) = parin;

for j=1:numel(parin)
  if isempty(parin{j})
    varout(j) = varin(j);
  end
end

for j = numel(parin)+1:numel(varin)
  varout(j) = varin(j);
end

if nargout == 1 & numel(varout) > 1
  varargout = {varout};
else
  varargout = varout;
end

  
  