function result=mahalanobis(x,locvct,varargin)

%MAHALANOBIS computes the distance of each observation in x
%   from the location estimate (locvct) of the data, 
%   relative to the shape of the data.  
%
% Required input arguments:
%                  x : data matrix (n observations in rows, p variables in columns)
%             locvct : location estimate of the data (p-dimensional vector)
%      cov or invcov : scatter estimate of the data or the inverse of the scatter estimate (pxp matrix)
%
% I/O: result=mahalanobis(x,locvct,'cov',covmat,'n',size(x,1),'p',size(x,2))
%   The user should only give the input arguments that have to change their default value.
%   The name of the input arguments needs to be followed by their value.
%   The order of the input arguments is of no importance.
%
% Examples:
%   result=mahalanobis(x,loc,'cov',covx,'n',10)
%   result=mahalanobis(x,loc,'p',2,'invcov',invcovx)
%   result=mahalanobis(x,loc,'invcov',invcovx)
%
% Output:
%   A vector containing the distances of all the observations to locvct.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Katrien Van Driessen
% Revisions by Sabine Verboven
% Last update on 18/09/2003
%

%Initialisation
n = size(x,1);
p = size(x,2);
if nargin<3
    error('Missing a required input variable')
end
counter=1;
default=struct('cov',0,'invcov',NaN);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
%reading the user's input 
if nargin>2
    %
    %placing inputfields in array of strings
    %
    for j=1:nargin-2
        if rem(j,2)~=0
            chklist{i}=varargin{j};
            i=i+1;
        end
    end 
    %
    %Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
    %
    while counter<=IN 
        index=strmatch(list(counter,:),chklist,'exact');
        if ~isempty(index) %in case of similarity
            for j=1:nargin-2 %searching the index of the accompanying field
                if rem(j,2)~=0 %fieldnames are placed on odd index
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
if size(locvct,2)==1
    locvct=locvct'; %converting to a rowvector 
end
if options.cov==0 & options.invcov==0
    error('The scatter matrix or its inverse is a required input argument.')
end
%%%%%%MAIN%%%%%%%%%
if ~isnan(options.invcov)
    covmat=options.invcov;
    if min(size(covmat))==1
        covmat=diag(covmat);
    end
else
    if min(size(options.cov))==1
        options.cov=diag(options.cov);
    end
    covmat=pinv(options.cov);
end
hlp=x-repmat(locvct,n,1); 
dist=sum(hlp*covmat.*hlp,2)'; 
result=dist;