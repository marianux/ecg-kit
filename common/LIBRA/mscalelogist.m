function result=mscalelogist(x,varargin)

%MSCALELOGIST is an M-estimator of scale with psi-function equal to
% (exp(x)-1)/(exp(x)+1) and with auxiliary location estimate.
% It is iteratively computed. It can resist 50% outliers. 
% If x is a matrix, the scale estimate is computed on the columns of x. The
% result is then a row vector. If x is a row or a column vector,
% the output is a scalar.
%
% The estimator is introduced in:
%   Rousseeuw, P.J. and Verboven, S. (2002),
%   "Robust estimation in very small samples", 
%   Computational Statistics and Data Analysis, 40, 741-758.
%
% Required input argument:
%    x: either a data matrix with n observations in rows, p variables in columns
%       or a vector of length n.
%
% Optional input arguments:
%      k: number of iteration steps (default value = 50)
%    sca: a starting value for the scale estimate
%         default value = 'madc'
%         other possibilities: 'qn'/'adm'/... 
%    loc: an auxiliary location estimate
%         default value = 'median'
%         other possibilities: 'hl'/'mlochuber'/...
%
% I/O: result=mscalelogist(x,'k',50,'sca','madc','loc','median')
%
% Example: result=mscalelogist(x,'sca','qn','k',10,'loc','hl')
%          result=mscalelogist(x,'loc','hl')
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
%
%Written by S.Verboven
%Revisions by Nele Smets
%Last update: 31/08/03 

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
    options.k=options.k;
    options.sca=options.sca;
    options.loc=options.loc;
end

if n==1 & p==1
    out=0;  %when X is a one by one matrix, all scale estimators must equal to 0
    return
elseif  n==1      
    x=x';   %we only want to work with column vectors
    n=p;
    p=1;
end

b=0.3739;
beta=0.5; %0.500038854875226   %quadl('(tanh(x./(2*b)).^2).*normpdf(x,0,1)',-5,5,1.e-15) ; denumenator of the scale estimator = constant
for i=1:p
    X=x(:,i); 
    j=1;
    s_0=feval(options.sca,X);
    t_0=feval(options.loc,X);
    step=s_0;
    if s_0==0
       out(i)=0;
    else
        while j<=options.k
            u=(X-t_0)/step;
            uu=tanh(u/(2*b)).^2;
            step=step*sqrt(sum(uu)/(n*beta));
            j=j+1;
        end
        out(i)=step;
        u=[];
    end
end

result = out;


