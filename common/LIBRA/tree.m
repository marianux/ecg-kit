function tree(objectorder,heights)

%TREE creates a tree in which the leaves represent
%   objects.  The vertical coordinate of the junction
%   of two branches is the dissimilarity between the
%   corresponding clusters (maximal 30 objects allowed).
%
% The algorithm is fully described in:
%   Kaufman, L. and Rousseeuw, P.J. (1990),
%   "Finding groups in data: An introduction to cluster analysis",
%   Wiley-Interscience: New York (Series in Applied Probability and
%   Statistics), ISBN 0-471-87876-6.
%
% Required input arguments:
%   objectorder :  order of objects
%   heights     : diameter of cluster before dividing it
%                 (=length of banner)
%
% I/O:
%   tree(objectorder,heights)
%
% Example (subtracted from the referenced book)
%   load agricul.mat
%   result = diana(agricul,[4 4],0,0,1);
%   tree(result.objectorder,result.heights)
%
% The output of TREE is a figure containing the
%   agglomerative (agnes) or divise (diana) tree.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Wai Yan Kong (May 2006)
% Last Revision: 28/09/2006


clf reset
whitebg([1 1 1]);

if (nargin<2)
    error('Two input arguments required')
elseif (nargin>2)
    error('Too many input arguments')
end

if(size(objectorder,2)~=size(heights,2)+1)
    error('Missing values in objectorder or heights')
end

Heights=H(heights);

number=size(objectorder,2);

if (number>30)
    error('Only 30 objects allowed')
end

%midden=[];
middle(1)=0;
Maxi=0;
Mini=0;
PrevMID=middle(1);
lengt=[]; %length
high=[]; %height
indices=[];
L=[];

[maxim,index]=max(Heights);
if(1<=index-1)
    [Prevmax,Previndex]=max(Heights(1:(index-1)));
else
    Prevmax=maxim;
    Previndex=index;
end
if(index+1<=number-1)
    [Postmax,Postindex]=max(Heights((index+1):(number-1)));
    Postindex=Postindex+index+1-1;
else
    Postmax=maxim;
    Postindex=index;
end

if(Postindex+1<=number-1)
    [Post2max,Post2index]=max(Heights(Postindex+1:number-1));
    Post2index=Post2index+Postindex+1-1;
else
    Post2max=Postmax;
    Post2index=Postindex;
end

if(Postindex+1<=Post2index-1)
    [Betweenmax,Betweenindex]=max(Heights(index+1:Postindex-1));
    Betweenindex=Betweenindex+index+1-1;
else
    Betweenmax=Postmax;
    Betweenindex=Postindex;
end

L=cat(2,L,maxim);

high=cat(2,high,[Post2max,Postmax,Betweenmax,maxim,Prevmax]);
indices=cat(2,indices,[Post2index,Postindex,Betweenindex,index,Previndex]);

M=0;
extra=number/2;
lengt(1)=Postindex-Previndex+extra+2;
PrevLEN=lengt(1);
rectangle('Position',[middle(1)-(lengt(1)/2),maxim,lengt(1),0.0001]);

NbanFirst=1;
k=1;
branch=0;
over=0;
Element=0;
Special=0;
Spec=0;
Sp=0;
S=0;
ww=0;
if(maxim<20)
    extrawaystick=0.8;
    extrawaytext=1.5;
else
    extrawaytext=0;
    extrawaystick=0;
end

right=0;
direct11=0;
direct111=0;
righttree=0;
if(index==1)
    direct1=1;
else
    direct1=0;
end

Last=size(lengt,2);
LastM=size(middle,2);
LastL=size(lengt,2);

while(k<=number)
    w=0;
    LastH=size(high,2);
    if((high(LastH)~=high(LastH-1)) & Element==0)
        while(branch==0)
            S=0;
            w=w+1;
            ww=ww+1;

            if(over==0)
                middle=cat(2,middle,middle(LastM)-(lengt(LastL)/2));
            end

            if(righttree==1)
                middle=cat(2,middle,middle(LastM)-(lengt(LastL)/2));
                righttree=0;
            end

            if(Sp==1)
                LastM=size(middle,2);
                LastL=size(lengt,2);
                middle(LastM)=middle(LastM)+lengt(LastL)/2;
                Sp=0;
                S=1;
            end
            LastM=size(middle,2);

            Output=maxim-Prevmax;
            rectangle('Position',[middle(LastM),Prevmax,0.0001,Output]);

            if(Betweenmax==Postmax | over==0)
                Post2max=Postmax;
                Post2index=Postindex;
            end

            if((Betweenmax==Postmax | over==0) & (right~=1 | over==0))
                Postmax=maxim;
                Postindex=index;
            end

            maxim=Prevmax;
            index=Previndex;

            if(NbanFirst<=(index-1))
                [Prevmax,Previndex]=max(Heights(NbanFirst:(index-1)));
                Previndex=Previndex+NbanFirst-1;
            else
                Prevmax=maxim;
                Previndex=index;
            end

            if(index+1<=Postindex-1)
                [Betweenmax,Betweenindex]=max(Heights(index+1:Postindex-1));
                Betweenindex=Betweenindex+index+1-1;
            else
                Betweenmax=Postmax;
                Betweenindex=Postindex;
            end

            if(right==1 & over~=0)
                high=cat(2,high,[Prevmax Prevmax]);
                indices=cat(2,indices,[Previndex Previndex]);
            else
                high=cat(2,high,Prevmax);
                indices=cat(2,indices,Previndex);
            end

            if(Postindex-Previndex<=0)
                lengt=cat(2,lengt,1);
            else
                if(Postindex-Previndex-direct111>0)
                    if(ww==1)
                        lengt=cat(2,lengt,Postindex-Previndex-direct111+2);
                    else
                        lengt=cat(2,lengt,Postindex-Previndex-direct111);
                    end
                else
                    if(ww==1)
                        lengt=cat(2,lengt,Postindex-Previndex+2);
                    else
                        lengt=cat(2,lengt,Postindex-Previndex);
                    end
                end
            end
            L=cat(2,L,maxim);

            if(direct111~=0)
                direct111=0;
            end

            LastH=size(high,2);
            LastL=size(lengt,2);
            LastM=size(middle,2);


            rectangle('Position',[middle(LastM)-(lengt(LastL)/2),maxim,lengt(LastL),0.0001]);

            if(over~=0)
                over=0;
            end

            if(right~=0)
                right=0;
            end

            if(NbanFirst==index)
                branch=1;
            end
        end

        rectangle('Position',[middle(LastM)-(lengt(LastL)/2),maxim-1+extrawaystick,0.0001,1-extrawaystick]);
        if(objectorder(k)>10)

            text(middle(LastM)-(lengt(LastL)/2)-0.25,maxim-2.0+extrawaytext,num2str(double(objectorder(k))));
            T=middle(LastM)-(lengt(LastL)/2-0.25);
        else

            text(middle(LastM)-(lengt(LastL)/2)-0.1,maxim-2.0+extrawaytext,num2str(double(objectorder(k))));
            T=middle(LastM)-(lengt(LastL)/2-0.1);
        end

        Maxi=max(Maxi,T);
        Mini=min(Mini,T);

        if(Betweenmax~=Postmax)
            high=high(1:LastH-2);
            indices=indices(1:LastH-2);

            Prevmax=Betweenmax;
            Previndex=Betweenindex;

            high=cat(2,high,Prevmax);
            indices=cat(2,indices,Previndex);
            LastH=size(high,2);

            LastL=size(lengt,2);
            middle(LastM)=middle(LastM)+lengt(LastL)/2;
            over=1;
        end
        direct1=0;
        Special=0;
        Spec=0;
    else
        if(Special==1)
            Special=0;
            Spec=1;
        end

        if(Sp==1)
            Sp=0;
        end

        if(S==1)
            S=2;
        end

        direct111=0;
        LastH=size(high,2);
        LastM=size(middle,2);
        LastL=size(lengt,2);

        if(direct1==0)

            rectangle('Position',[middle(LastM)+(lengt(LastL)/2),high(LastH)-1+extrawaystick,0.0001,1-extrawaystick]);

            text(middle(LastM)+lengt(LastL)/2-0.1,high(LastH)-2.0+extrawaytext,num2str(double(objectorder(k))));

            T=middle(LastM)+lengt(LastL)/2-0.1;
            Maxi=max(Maxi,T);
            Mini=min(Mini,T);
        else

            rectangle('Position',[middle(LastM)-(lengt(LastL)/2),high(LastH)-1+extrawaystick,0.0001,1-extrawaystick]);

            text(middle(LastM)-lengt(LastL)/2-0.1,high(LastH)-2.0+extrawaytext,num2str(double(objectorder(k))));

            T=middle(LastM)-lengt(LastL)/2-0.1;
            Maxi=max(Maxi,T);
            Mini=min(Mini,T);
            direct1=0;
            direct11=1;
        end

        if(high(LastH-1)==high(LastH))
            high=high(1:LastH-2);
            indices=indices(1:LastH-2);
        else
            high=high(1:LastH-1);
            indices=indices(1:LastH-1);
        end
        LastH=size(high,2);

        maxim=high(LastH);
        index=indices(LastH);

        if(LastH-1>=1)
            Postmax=high(LastH-1);
            Postindex=indices(LastH-1);
        end

        if(LastH-2>=1)
            Post2max=high(LastH-2);
            Post2index=indices(LastH-2);
        end

        if(index+1<=Postindex-1)
            [Betweenmax,Betweenindex]=max(Heights(index+1:Postindex-1));
            Betweenindex=Betweenindex+index+1-1;
        else
            Betweenmax=Postmax;
            Betweenindex=Postindex;
        end

        if(Betweenmax~=Postmax)
            Prevmax=Betweenmax;
            Previndex=Betweenindex;
        else
            Prevmax=maxim;
            Previndex=index;
        end

        if(high(LastH)>high(LastH-1))
            righttree=1;
        end

        high(LastH)=Prevmax;
        indices(LastH)=Previndex;

        LastM=size(middle,2);
        LastL=size(lengt,2);
        Last=size(L,2);

        if(LastM-1>=1)
            middle=middle(1:LastM-1);
        end

        if(LastL-1>=1)
            lengt=lengt(1:LastL-1);
        end

        if(Last-1>=1)
            L=L(1:Last-1);
        end

        LastM=size(middle,2);
        LastL=size(lengt,2);
        Last=size(L,2);

        if(Last>0)
            while(maxim>L(Last))
                lengt=lengt(1:Last-1);
                L=L(1:Last-1);
                Last=size(L,2);
            end
        end
        LastL=size(lengt,2);

        middle(LastM)=middle(LastM)+(lengt(LastL)/2);

        if(Betweenmax~=Postmax)
            Element=0;
        else
            Element=1;
            middle(LastM)=middle(LastM)-(lengt(LastL)/2);
        end

        if(righttree==1 | direct11==1)
            if(Spec==1)
                Sp=1;
            end

            PrevMAX=L(1);
            [maxim,index]=max(Heights(NbanFirst+1:number-1));
            index=index+NbanFirst+1-1;

            if(NbanFirst+1<=index-1)
                [Prevmax,Previndex]=max(Heights(NbanFirst+1:index-1));
                Previndex=Previndex+NbanFirst+1-1;
            else
                Prevmax=maxim;
                Previndex=index;
            end

            if(index+1<=number-1)
                [Postmax,Postindex]=max(Heights(index+1:number-1));
                Postindex=Postindex+index+1-1;
            else
                Postmax=maxim;
                Postindex=index;
            end

            if(index+1<=Postindex-1)
                [Betweenmax,Betweenindex]=max(Heights(index+1:Postindex-1));
                Betweenindex=Betweenindex+index+1-1;
            else
                Betweenmax=Postmax;
                Betweenindex=Postindex;
            end

            if(Postindex+1<=number-1)
                [Post2max,Post2index]=max(Heights(Postindex+1:number-1));
                Post2index=Post2index+Postindex+1-1;
            else
                Post2max=Postmax;
                Post2index=Postindex;
            end


            rectangle('Position',[PrevMID+(PrevLEN/2),maxim,0.0001,PrevMAX-maxim]);

            LastM=size(middle,2);
            LastMID=PrevMID+(PrevLEN/2);
            PrevMID=LastMID;
            middle=[];
            middle(1)=LastMID;
            LastM=size(middle,2);

            high=[];
            indices=[];
            high=cat(2,high,[Post2max,Postmax,Betweenmax,maxim,Prevmax]);

            indices=cat(2,indices,[Post2index,Postindex,Betweenindex,index,Previndex]);
            LastH=size(high,2);

            L=[];
            L=cat(2,L,maxim);
            Last=size(L,2);
            lengt=[];
            if(Postindex>Previndex)
                if(M<extra)
                    M=M+1;
                end
                lengt(1)=Postindex-Previndex+(extra-M);

                rectangle('Position',[LastMID-lengt(1)/2,maxim,lengt(1),0.0001]);
                direct1=1;
            elseif(Postindex<=Previndex & Betweenmax==Postmax)
                lengt(1)=1;
                rectangle('Position',[LastMID-1/2,maxim,1,0.0001]);
                %%%%%%%

                rectangle('Position',[LastMID-1/2,maxim-1+extrawaystick,0.0001,1-extrawaystick]);
                if(objectorder(k+1)>10)

                    text(LastMID-1/2-0.25,maxim-2.0+extrawaytext,num2str(double(objectorder(k+1))));
                    Mini=min(Mini,LastMID-1/2-0.25);
                else

                    text(LastMID-1/2-0.1,maxim-2.0+extrawaytext,num2str(double(objectorder(k+1))));
                    Mini=min(Mini,LastMID-1/2-0.1);
                end

                rectangle('Position',[LastMID+1/2,maxim-1+extrawaystick,0.0001,1-extrawaystick]);

                text(LastMID+1/2-0.1,maxim-2.0+extrawaytext,num2str(double(objectorder(k+2))));
                Maxi=max(Maxi,LastMID+1/2-0.1);

                k=k+2;
            end
            PrevLEN=lengt(1);
            LastL=size(lengt,2);

            if(Element==1)
                Element=0;
            end

            if(direct11==1)
                if(NbanFirst+1<index)
                    middle=cat(2,middle,middle(LastM)-lengt(LastL)/2);
                    LastM=size(middle,2);
                    direct111=1;
                    over=1;
                end
                direct11=0;
            end

            Special=1;
            PrevMAX=maxim;
        end
        if(S==2 && k+1==number)

            rectangle('Position',[LastMID+lengt(1)/2,PrevMAX-1+extrawaystick,0.0001,1-extrawaystick]);

            text(LastMID+lengt(1)/2-0.1,PrevMAX-2.0+extrawaytext,num2str(double(objectorder(k+1))));
            Maxi=max(Maxi,LastMID+lengt(1)/2-0.1);
            Mini=min(Mini,LastMID+lengt(1)/2-0.1);
            k=k+1;
        end
        over=1;
    end
    branch=0;


    k=k+1;
    NbanFirst=NbanFirst+1;


end
axis([Mini-1,Maxi+1,min(Heights)-5,max(Heights)+1]);
XT=[];
set(gca,'XTick',XT);
set(gca,'XTickLabel',[]);

%---
function res=H(vector)

lengt=size(vector,2);
Nvector=vector;
i=1;
for i=1:lengt
    for j=i+1:lengt
        if (vector(i)==vector(j))
            Nvector(j)=vector(j)+i*0.0001;
            i=i+1;
        end
    end
end
res=Nvector;


