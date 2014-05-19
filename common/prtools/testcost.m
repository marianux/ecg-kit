function e = testcost(x,w,C,lablist)
%TESTCOST compute the error using the cost matrix C
%
% Be Aware! This routine needs a lot of maintenance!!   (RD)
%
%   E = TESTCOST(A,W,C,LABLIST)
%   E = TESTCOST(A*W,C,LABLIST)
%   E = A*W*TESTCOST([],C,LABLIST)
%
% INPUT
%   A       Dataset
%   W       Trained classifier mapping
%   C       Cost matrix
%   LABLIST Labels corresponding to the entries of C
%
% OUTPUT
%   E       Total classification cost
%
% DESCRIPTION
% Compute the misclassification cost using the cost matrix C. In
% LABLIST the corresponding classes should be listed. Note that
% classifier W should be a 'genuine' classifier, in the sense that the
% output of the classifier should be transformed using CLASSC. See
% for an example the help of costm.m.
%
% SEE ALSO
% testd, costm, setcost, getcost

% Copyright: D.M.J. Tax, davidt@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

if nargin < 4
	lablist = [];
end
if nargin < 3
	C = [];
end

% Do it the same as testc:
% When no input arguments are given, we just return an empty mapping:
if (nargin == 0) | (isempty(x))
	
	% No mercy: all parameters have to be defined:
	if isempty(C)
		error('Lablist should be defined for empty mapping');
	end
	if isempty(w)
		error('C-matrix should be defined for empty mapping');
	end
		
	% Sometimes Prtools is crazy, but fun!:
	e = prmapping(mfilename,'fixed',{w,C});
	return

elseif (~isempty(w)) & (ismapping(w))

	% we first have to map the data using the classifier:
	ismapping(w);
	istrained(w);

	e = feval(mfilename,x*w,C,lablist);

else
	% Now we are doing the actual work:
   % find out where the cost-matrix and lablist are:
	lablist = C;
	C = w;
   if isempty(C)
      cost = getuserdata(x);
      C = cost.C;
      lablist = cost.lablist;
   end
	% Check if everything is defined:
	if isempty(C)
		error('Cost matrix is not defined!');
	end
	if isempty(lablist)
		error('Lablist is not defined!');
	end
	% true labels:
	[nlab1,lablist1]=getnlab(x);
	% match with C-lablist
	I1 = matchlablist(lablist1,lablist);
	% give a warning if it does not fit:
	if any(I1==0)
		error('Some objects have labels which are not defined in C!');
	end
	% and fix it:
	nlab1 = I1(nlab1);

	% estimated labels:
	[mx,nlab2] = max(x,[],2);
	lablist2 = getfeatlab(x);
%La= lablist2(nlab2,:);
	% match with C-lablist
	I2 = matchlablist(lablist2,lablist);
	% give a warning if it does not fit:
	if any(I2==0)
		error('The classifier outputs a label which are not defined in C!');
	end
	% and fix it:
	nlab2 = I2(nlab2);
%Lb= lablist(nlab2,:);
%[La Lb]

	% finally, we can compute the cost:
	e = 0;
	n = length(nlab1);
	for i=1:n
		%DXD: to make the cost matrix definition consistent with that one
		%used in costm.m, this should be the order of the lab1 and lab2:
		e = e + C(nlab1(i),nlab2(i));
	end

	% return the average cost (i.e. normalized by the number of objects):
	e = e/n;

end

return
