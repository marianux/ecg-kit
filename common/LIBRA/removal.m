function restX=removal(X,r,k)

%REMOVAL deletes rows(r) or columns(k) from X
% Type in zero if you do not want to delete any rows (or columns)
% r= rowvector of object numbers you want to remove 
% k= rowvector of variable numbers you want to remove
%
%I/O: restX=removal(X,r,k)
%
%Created by S.Verboven (1999)
%Last edited 25/04/2000

if nargin <3
   error('All inputarguments must be given. Type in zero if you do not want to delete any rows(or columns)')
end
restX=X;
if size(r,1)~=1
   r=r';
end
if size(k,1)~=1
   k=k';
end
if k==0 & r~=0
   restX(r,:)=[];
end
if r==0 & k~=0
   restX(:,k)=[];
end
if r~=0 & k~=0
   restX(r,:)=[];
   restX(:,k)=[];
end
