%PRLOAD Similar to LOAD matfile, but converts old datasets, datafiles, mappings
%
%  PRLOAD(FILE)
%  S = PRLOAD(FILE)

function out = prload(file)

warning('off','MATLAB:unknownObjectNowStruct');
warning('off','MATLAB:indeterminateFields');
warning('off','MATLAB:unknownElementsNowStruc');
warning('off','MATLAB:elementsNowStruc');

ppath = path;
ss = which(['data' 'set']); % avoid conversion to prdataset in prtools5
if ~isempty(ss) & isempty(strfind(ss,'prtools'))
  rmpath(fileparts(fileparts(ss)));
end

s = load(file);

path(ppath);

if ~isstruct(s)
  return
end  

fields = fieldnames(s);
for j=1:numel(fields)
  x = getfield(s,fields{j});
  if isstruct(x)
    if isfield(x,'mapping_file')
      x = prmapping(x);
    elseif isfield(x,'rootpath')
      x = prdatafile(x);
    elseif isfield(x,'labtype')
      x = prdataset(x);
    elseif isfield(x,'ll')
      x = prdataset(x);
    end
  end
  s = setfield(s,fields{j},x);
end

if nargout > 0
  out = s;
else
  for j=1:numel(fields)
    assignin('caller',fields{j},getfield(s,fields{j}));
  end
end
  

