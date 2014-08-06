%CLASSC Convert classifier to normalized classifier (yielding confidences)
%
%  V = CLASSC(W)
%  V = W*CLASSC
%  D = CLASSC(A*W)
%  D = A*W*CLASSC
%  D = CLASSC(A,W)
%
% INPUT
%  W  Trained or untrained classifier
%  A  Dataset
%	
% OUTPUT
%  V  Normalized classifier producing confidences instead of 
%     densities or distances (after training if W is untrained)
%
% DESCRIPTION
% The trained or untrained classifier W may yield densities or unnormalised
% confidences. The latter holds for two-class discriminants like FISHERC
% and SVC as well as for neural networks. Such classifiers use or should
% use CNORMC to convert distances to confidences. In multi-class problems
% as well as in combining schemes they do not produce normalises
% confidences. These outcomes, like the density outcomes of classifiers
% liek QDC, LDC and PARZENC, can be converted by CLASSC into confidences:
% the sum of the outcomes will be one for every object.
%
% In case W is a one-dimensional mapping, it is converted into a two-class
% classifier, provided that during the construction a class label was 
% supplied. If not, the mapping cannot be converted and an error is
% generated.
%
% CLASSC lists the outcomes on the screen in case no output argument is
% supplied. Also true and estimated labels are supplied.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, CNORMC, LABELD

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: classc.m,v 1.5 2010/02/23 15:21:54 duin Exp $

function w = classc(w,flag)

		if nargin < 2, flag = 0; end % flag forces non=combiner behavior avoiding recursion
	
	if (nargin == 0)

		% Untrained mapping.
		w = prmapping('classc','combiner',flag);			

	elseif (ismapping(w))

		% If mapping is stacked or parallel, recurse over the individual
		% sub-mappings and call CLASSC for each of them.
		if ((isstacked(w)) | (isparallel(w))) & (flag == 0)
			v = cell(1,length(w.data));
			for j = 1:length(w.data)
				if ismapping(w.data{j}) % the parallel combiner may have nonmapping data
					v{j} = feval(mfilename,w.data{j});
				else
					v{j} = w.data{j};
				end
			end
			w = setdata(w,v);
			w = feval(mfilename,w,1);  % and here CLASSC is called for the combiner avoiding recursion
		else
			conv = get(w,'out_conv');
			if (conv < 1)
				% Set the "normalization" bit in the mapping's output conversion flag
				w = set(w,'out_conv',conv+2);			
			end;
		end
	elseif (isdataset(w))
		if ismapping(flag)
			if nargout == 1
				w = feval(mfilename,w*flag);
			else
				feval(mfilename,w*flag);
				clear w;
			end
			return
		end
		w = w*normm;
		w = w*costm;
		if nargout == 0 % list outcomes on the screen
			ww = +w;
			ss = repmat('-',1,9*size(ww,2));
			fprintf('\n True   Estimated        Class \nLabels    Labels      Confidences\n');
			fprintf('------------------%s\n',ss);
			nlab = getnlab(w);
			[wmax,K] = max(ww,[],2);
			lablist = getlablist(w);
			if ~isempty(lablist) & ~ischar(lablist)
				nlab = lablist(nlab);
				K = lablist(K);
			end
			for j=1:size(ww,1)
				if (nlab(j) ~= K(j))
					fprintf(' %3.0f    ->%3.0f    ',nlab(j),K(j));
				else
					fprintf(' %3.0f      %3.0f    ',nlab(j),K(j));
				end
				fprintf('  %7.4f',ww(j,:));
				fprintf('\n');
			end
			lablist = getlablist(w);
			if ischar(lablist)
				fprintf('\n');
				for j=1:size(lablist,1)
					fprintf('   %i %s\n',j,lablist(j,:));
				end
			end
				
				
			clear w;
		end
	else
		error('input should be mapping or dataset');
	end

return
