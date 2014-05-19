%READDATAFILE Read one of the datafiles
%
%    [B,NEXT,J] = READDATAFILE(A,N)
%
% INPUT     
%   A           Datafile
%   N           Number of the file to be read
%
% OUTPUT
%   B           Dataset stored in file N
%   NEXT        Number of next file to be read, 0 if done
%   J           Indices of objects in A
%
% DESCRIPTION
% A datafile points to a dataset stored in a series of files. This
% routine reads one of them, but is designed to read them all in a loop.
% A typical example is shown below, computing the overall mean per class.
% If the preprocessing field of A is set, the listed preprocessing is
% applied before returning.
% If the mappings field of A is set, the listed mappings are applied
% to B before returning.
%
% As the objects in A may be randomly distributed over the files, a 
% reordering is performed internally in this routine. Consequently,
% objects may be returned in a different order than stored in A.
%
% [m,k,c] = getsize(a);
% nobjects = classsizes(a);
% u = zeros(c,k);
% next = 1;
% while next > 0
%    [b,next] = readdatafile(a,next)
%    u = u + meancov(b) .* repmat(nobjects',1,k);
% end
% u = u ./ repmat(classsizes(a)',1,k);
%
% SEE ALSO
% PRDATASET, DATAFILE
