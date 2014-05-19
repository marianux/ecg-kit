function result = robstd(x,varargin)

%ROBSTD standardizes a matrix X columnwise by substracting a robust estimate
% of the location and by dividing through a robust estimate of the scale.
%
% Required input arguments:
%          x: Data matrix. Rows of x represent observations, and columns represent 
%             variables.
%
% Optional input arguments:
%   loc     : Name of the m-file which computes the location estimate (default='median')
%   sca     : Name of the m-file which computes the scale estimate (default='madc')
%   iterloc : If one uses an iterated M-estimator of location, iterloc contains the
%             number of iteration steps.
%   itersca : If one uses an iterated M-estimator of scale, iterloc contains the
%             number of iteration steps.
%         h : If one uses the unimcd estimator; h is the quantile of observations 
%             whose variance will be minimized.  
%             Any value between n/2 and n may be specified. The default value is [0.75*n].
%
% I/O: result=robstd(x,'loc','unimcd','h',20)
%
% The ouput is a structure containing the following fields:
%       loc : vector with the location estimates of each variable.
%       sca : vector with the scale estimates of each variable.
%      xstd : a matrix like x, containing the standardized observations
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Nele Smets on 25/03/2004
% Last update: 08/04/2004

[n,p]=size(x);
%
% initialization with defaults
%
alfa=0.75;
hdef=min(floor(2*floor((n+p+1)/2)-n+2*(n-floor((n+p+1)/2))*alfa),n);
counter=1;
default=struct('sca','madc','loc','median','iterloc',50,'itersca',50,'h',hdef);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
%
if nargin>1
    %
    % placing inputfields in array of strings
    %
    for j=1:nargin-1  
        if rem(j,2)~=0
            chklist{i}=varargin{j};    
            i=i+1;
        end
    end 
    %
    % Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
    %
    while counter<=IN 
        index=strmatch(list(counter,:),chklist,'exact');
        if ~isempty(index) % in case of similarity
            for j=1:nargin-1 % searching the index of the accompanying field
                if rem(j,2)~=0 % fieldnames are placed on odd index
                    if strcmp(chklist{index},varargin{j})
                        I=j;
                    end
                end
            end
            options=setfield(options,chklist{index},varargin{I+1});
            index=[];
        end
        counter=counter+1;
    end
end
sca=options.sca;
loc=options.loc;
iterloc=options.iterloc;
itersca=options.itersca;
h=options.h;

if ~ismember(sca,{'qnm','madc','adm','mscalelogist','unimcd'})
    error('No appropriate scale estimator was given.')
end
if ~ismember(loc,{'median','hl','mloclogist','mlochuber','unimcd'})
    error('No appropriate location estimator was given.')
end

if strmatch(loc,'unimcd')&strmatch(sca,'unimcd')
    [out.loc,out.sca]=unimcd(x,h);
else
    switch loc
    case {'hl','median'} 
        out.loc=feval(loc,x);
    case 'unimcd'
        out.loc=feval(loc,x,h);
    case 'mlochuber' 
        out.loc=feval(loc,x,'k',iterloc);
    case 'mloclogist'
        out.loc=feval(loc,x,'k',iterloc);
    end
    switch sca
    case {'adm', 'qnm','madc'} 
        out.sca=feval(sca,x);
    case 'unimcd'
        [l,out.sca]=feval(sca,x,h);
    case 'mscalelogist'
        out.sca=feval(sca,x,'k',itersca);
    end
end
out.xstd=(x-repmat(out.loc,n,1))./repmat(out.sca,n,1);

result = out;