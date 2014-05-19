function result=mlochuber(x,varargin)

%MLOCHUBER is an M-estimator of location with psi-function equal to
% min(1.5,max(-1.5,x)) and with an auxiliary scale estimate.
% It is iteratively computed. It can resist 50% outliers. 
% If x is a matrix, the location estimate is computed on the columns of x. The
% result is then a row vector. If x is a row or a column vector,
% the output is a scalar.
%
% The estimator is described in:
%   Huber, P. (1981), Robust Statistics, Wiley, New York. 
% Its behavior in small samples is discussed in:
%   Rousseeuw, P.J. and Verboven, S. (2002),
%   "Robust estimation in very small samples", 
%   Computational Statistics and Data Analysis, 40, 741-758.
%
% Required input argument:
%    x: either a data matrix with n observations in rows, p variables in columns
%       or a vector of length n.
%
% Optional input arguments:
%    k: number of iteration steps (default value = 50)
%  loc: a starting value for the location estimate
%       default value = 'median'
%       other possibilities : 'hl'/'mloclogist'/...
%  sca: an auxiliary scale estimate
%       default value = 'madc'
%       other possibilities: 'qn'/'adm'/'mscalelogist'/...
%
% I/O: result=mlochuber(x,'k',50,'loc','median','sca','mad')    
%
% Examples: result=mlochuber(x,'loc','mlochuber','sca','qn');
%           result=mlochuber(x,'sca','mscalelogist','k',10);
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by S. Verboven
% Last update 12/02/04 

[n,p]=size(x);
%
% initialization with defaults
%
counter=1;
sca='madc';
loc='median';
default=struct('k',50,'sca',sca,'loc',loc);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
%
%
if nargin>1
    %
    % placing inputfields in array of strings
    %
    for j=1:nargin-1  
        if rem(j,2)~=0
            chklist{i}= varargin{j};    
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
    options.k=options.k;
    options.sca=options.sca;
    options.loc=options.loc;  
end

if n==1 & p==1
    out=x;     %when X is a one by one matrix, all location estimators must be equal to that matrix
    return
elseif n==1
    x=x';      %we only want to work with column vectors
    n=p;
    p=1;
end

if n==2  % all location estimators must equal the average for n=2
    out=mean(x,1);
    return
end

alfa=0.866385597462284; %2*normcdf(1.5,0,1)-1; constant denumenator
for i=1:p  
    X=x(:,i);
    t_0=feval(options.loc,X);
    tstep=t_0;
    s_0=feval(options.sca,X);
    if (s_0==0)
        out(i)=t_0;
    else
         j=1;
         while j<=options.k  
            z=(X-tstep)/s_0;
            y(abs(z)<=1.5)=z(abs(z)<=1.5);
            y(abs(z)>1.5)=1.5*sign(z(abs(z)>1.5));
            tstep=tstep+s_0*(sum(y)/(n*alfa)); %updating location estimate
            y=[];  %clear helpvector
            j=j+1;
         end
         out(i)=tstep;
    end
end

result = out;
 

