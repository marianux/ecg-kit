function result=mona(x,banner)

%MONA returns clusters according to mona (monothetic
%clustering)
%
%The algorithm is fully described in:
%   Kaufman, L. and Rousseeuw, P.J. (1990),
%   "Finding groups in data: An introduction to cluster analysis",
%   Wiley-Interscience: New York (Series in Applied Probability and
%   Statistics), ISBN 0-471-87876-6.
%
% Required input argument:
%  x : Data matrix (rows = observations, columns = variables)
%      containing only binary values.
%      The missing values are indicated by the value 2
%
% Optional input argument:
%  banner : draws picture
%       Possible values are 0 : do not create a banner
%                           1 : create a banner
%       (in case banner is not given, banner = 0)
% I/O:
%   result=mona(x,1)
%
% Example (subtracted from the referenced book)
%   load animal.mat
%   result=mona(animal,1);
%
% The output of MONA is a structure containing
%   result.matrix     : revised matrix xx (contains only 0 and 1 and missing
%                       values are estimated)
%   result.number     : number of observations
%   result.var        : number of variables
%   result.ner        : order of objects
%   result.lava       : variable used for separation
%   result.nban       : seperation step
%                       (The value on the ith index is the stepnumber from
%                       the seperation of the elements with index 1 to i
%                       from the elements with index i+1 to result.number)
%                       If it equals zero, then there was no separation
%
% And MONA will create the plot banner if banner equals 1.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Wai Yan Kong (May 2006)

resl=[];

% Check whether the number of input is correct

if (nargin<1)
    error('One input argument required (datamatrix)')
elseif (nargin<2)
    banner=0;
elseif (nargin>2)
    error('Too many input arguments')
end

% Define number of observations and variables

S=size(x);
S1=S(1);
S2=S(2);

resl.number=S(1);
resl.var=S(2);

% Check whether x is a matrix containing only binary values
% and whether x has missing values
% We make a revised matrix xx

missing=zeros(1,S2);
for j = 1:S2
    for i= 1:S1
        if (x(i,j)==0 | x(i,j)==1)
            missing(j)=missing(j);
            xx(i,j)=x(i,j);
        elseif (x(i,j)==2)
            missing(j)=missing(j)+1;
            for k=1:S2
                a=0;
                b=0;
                c=0;
                d=0;

                if k~=j
                    if x(i,k) ~= 2

                        for t=1:S1
                            if (x(t,j)==1 & x(t,k)==1)
                                a=a+1;
                            elseif (x(t,j)==1 & x(t,k)==0)
                                b=b+1;
                            elseif (x(t,j)==0 & x(t,k)==1)
                                c=c+1;
                            elseif (x(t,j)==0 & x(t,k)==0)
                                d=d+1;
                            end
                        end
                    end
                end
                association(k)=abs(a*d-b*c);
                resinbetween(k)=a*d-b*c;
                [C,I]=max(association);

                if resinbetween(I)>0
                    xx(i,j)=x(i,I);
                elseif resinbetween(I)<0
                    xx(i,j)=1-x(i,I);
                end
            end
        else
            error('inputmatrix must have binary values or value 2 for missing values')
        end
    end
end

TotalMissing=sum(missing);
fprintf(1,'This inputmatrix has %d missing values\n',TotalMissing)

% Check situations where Mona is not applicable

One=0;
for j=1:S2
    if missing(j)>=1
        One=One+1;
    end
end
if One==S2
    error('each variable has at least one missing value')
end


for j=1:S2

    if (missing(j)>= (S1/2))
        fprintf(1,'Variable %d has %d missing values\n',j,missing(j))
        error('The number of missing values for some variable equals or is more than half of the number of objects')
    end


    gelijk=0;
    for s=1:S1-1
        if xx(s,j)==xx(s+1,j)
            gelijk=gelijk+1;
        end
    end
    if gelijk==S1-1
        fprintf(1,'Variable %d has identical values\n',j)
        error('all values are identical for some variable')
    end

end

% We make the rowvector kx by reading the revised matrix xx
% row by row

Lengte=S1*S2;
kx=zeros(1,Lengte);

for i = 1:S1
    for j= 1:S2
        kx((i-1)*S2+j) = xx(i,j);
    end
end


% Actual calculations

[ner,lava,nban]=monac(resl.number,resl.var,kx);

% We want lava and nban to be vectors of length S1-1

lava2=zeros(1,S1-1);
nban2=zeros(1,S1-1);

for i = 1:S1-1
    lava2(i) = lava(i+1);
    nban2(i) = nban(i+1);
end

% Create a banner

if (banner==1)
    whitebg([1 1 1]);
    stepmax=max(nban2)+1;
    Y=nban2;
    Y(find(Y==0))=stepmax;
    Y1=fliplr(Y);
    barh(Y1,1)
    title 'Banner of Mona'
    xlabel('Separation Step')
    XT=0:stepmax;
    set(gca,'XTick',XT);
    set(gca,'XTickLabel',XT);
    axis([0,stepmax,0.5,S1-0.5])
    ylabel('Objects')
    La1=num2str(lava2);
    La2=fliplr(La1);
    for k=1:S1-1
        La3(k)=La2(k+((k-1)*2));
    end
    for i=1:S1-1
        if(Y1(i)<stepmax)
            text(Y1(i)+0.3,i,La3(i))
        end
    end
    YT=0.5:S1;
    set(gca,'YTick',YT);
    set(gca,'YTickLabel',fliplr(ner));

elseif ((banner~=0) & (banner~=1) & (nargin==2))
    error('banner must equals 0 or 1')
end


% Putting things together

result = struct('Matrix',xx,'Observations',S1,'Variables',S2,...
    'Objectorder',ner,'UsedVariable',lava2,'SeperationStep',nban2);




