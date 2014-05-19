function [indexes,RM,segind]=signaltest(x,step,mode,th1)
%
% detecting severe quality loss in ECG
%   In:
%       x = evaluation signal (square of the WT in scale 1)
%       step = samples to consider in each subinterval for local evaluation (use samplig frequency, by default 1000)
%       mode = signaltest mode 1:RMS max test (default)
%
%   Out:
%   indexes - samples with possible quality loss
%   
%   RM -  root of the mean of signal x values
%   segind - samples corresponding to the RM evaluation (beginning and end of each segment, beginning of each subsegment)
%
%   The root of the mean of signal x (RM) in each 2^16 maximum length segment
%   is compared with the global segment threshould th1=2*step
%   If RM>th1 the RM in each subsegment of length step is compared with 
%   the the local threshould th2=step and indexes for which RM>th2 are retrieved   
%
%   Created by Rute Almeida (rbalmeid@unizar.es) on May 14, 2009.
%   Revision 10MAR2009


if nargin<4
    %th1=2*step*10;
    %th1=1.75*step*10;
    %th2=step*10;
    th1=sqrt(mean(x.^2));
    if nargin<3
        mode=1; 
        if nargin<2
            step=1000;
        end
    end
end
if isempty(mode)
end
seglength=length(x);

segment=[-seglength+1 0];
RM=[];detector=[];
X=zeros(size(x));
segind=0;
lead=1; %#ok<NASGU>
flag=0;
n=length(x);
while segment(end)<n
    segment=[(segment(1)+seglength) (min(segment(end)+seglength,length(x)))];
    X((segment(1)):(segment(end)))=x;
    count=ceil(length(x)/step);
    intervals=[1 fix(n/4); fix(n/4) fix(n/2); fix(n/2) fix(n*3/4); fix(n*3/4) length(x)];
    RM_seg=[sqrt(mean(x(intervals(1,1):intervals(1,2)).^2)) sqrt(mean(x(intervals(2,1):intervals(2,2)).^2)) sqrt(mean(x(intervals(3,1):intervals(3,2)).^2)) sqrt(mean(x(intervals(4,1):intervals(4,2)).^2))];
    RM_ind=RM_seg>(th1*1.2);
    th2=mean(RM_seg(~RM_ind));    
    if sum(RM_ind)>0
        flag=1;
        for j=1:count
            segind=[segind segment(1)+(j-1)*step]; %#ok<AGROW>
            RM=[RM sqrt(mean(x((segind(end)-segment(1)+1):min(segind(end)-segment(1)+1+step,length(x))).^2))]; %#ok<AGROW>
        end
     segind(1)=[];
    else
        segind= intervals;
        RM=[RM  RM_seg]; %#ok<AGROW>
    end
end
if flag==1;
detector=segind(find(RM>(2*th2))); %#ok<FNDSB>
end
if ~isempty(detector)
    indexes=ones(step,1)*detector+(ones(length(detector),1)*(1:step))';
    %indexes=(ones(step,1)*detector*(step-1)+(ones(length(detector),1)*(1:step))');
    indexes=indexes(:);
else
    indexes=[];
end

if length(indexes)>0.9*length(x(:,1))
    indexes=[];
end
