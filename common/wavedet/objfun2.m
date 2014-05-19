function f = objfun2(X,C,weight)
% error matrix to minimize: to find direction and point
%
% Rute Almeida 2.Dec.2004
% Last update:
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13

%X=[a ;v]
if nargin<3
    weight=ones(size(C));
end
a=NaN*ones(1,3);
v=NaN*ones(1,3);
a(1:size(X,2))=X(1,:);
v(1:size(X,2))=X(2,:);
%aux=[C-[C(:,1) v(2)/v(1)*(C(:,1)-a(1))+a(2)  v(3)/v(1)*(C(:,1)-a(1))+a(3)]];
aux=[C(:,1) v(2)/v(1)*(C(:,1)-a(1))+a(2)  v(3)/v(1)*(C(:,1)-a(1))+a(3)];
aux=weight.*(C-aux(:,1:size(C,2)));
f=aux(:);
%f =sqrt(sqrt(trace(aux'*aux)));