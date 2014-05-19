%CONFMAT Construct confusion matrix
% 
%  [C,NE,LABLIST] = CONFMAT(LAB1,LAB2,METHOD,FID)
%
% INPUT
%  LAB1        Set of labels
%  LAB2        Set of labels
%  METHOD      'count' (default) to count number of co-occurences in
%	             LAB1 and LAB2, 'disagreement' to count relative
%		           non-co-occurrence.
%  FID         Write text result to file
%
% OUTPUT
%  C           Confusion matrix
%  NE          Total number of errors (empty labels are neglected)
%  LABLIST     Unique labels in LAB1 and LAB2
%
% DESCRIPTION
% Constructs a confusion matrix C between two sets of labels LAB1 
% (corresponding to the rows in C) and LAB2 (the columns in C). The order of 
% the rows and columns is returned in LABLIST. NE is the total number of 
% errors (sum of non-diagonal elements in C).
%
% When METHOD = 'count' (default), co-occurences in LAB1 and LAB2 are counted 
% and returned in C. When METHOD = 'disagreement', the relative disagreement 
% is returned in NE, and is split over all combinations of labels in C
% (such that the rows sum to 1). (The total disagreement for a class equals
% one minus the sensitivity for that class as computed by TESTC).
%
%   [C,NE,LABLIST] = CONFMAT(D,METHOD)
%
% If D is a classification result D = A*W, the labels LAB1 and LAB2 are 
% internally retrieved by CONFMAT before computing the confusion matrix.
%
% When no output argument is specified, or when FID is given, the
% confusion matrix is displayed or written a to a text file. It is assumed
% that LAB1 contains true labels and LAB2 stores estimated labels.
%
% EXAMPLE
% Typical use of CONFMAT is the comparison of true and and estimated labels
% of a testset A by application to a trained classifier W: 
% LAB1 = GETLABELS(A); LAB2 = A*W*LABELD.
% More examples can be found in PREX_CONFMAT, PREX_MATCHLAB.
% 
% SEE ALSO
% MAPPINGS, DATASETS, GETLABELS, LABELD

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: confmat.m,v 1.7 2008/10/14 21:34:32 duin Exp $

function [error_idx, error_type_idx, error_type_lablist] = classification_errors(arg1)

	prtrace(mfilename);

	% Check arguments.
    if isdataset(arg1)
        lab1 = getlabels(arg1); 
        lab2 = arg1*labeld;
    else
        error('A result dataset should be provided.');
    end
		
	% Renumber LAB1 and LAB2 and find number of unique labels.

	m = size(lab1,1);
	if (m~=size(lab2,1))
		error('LAB1 and LAB2 have to have the same lengths.');
	end
	[nlab1,lablist1] = renumlab(lab1);
	[nlab2,lablist2] = renumlab(lab2);

%     n1 = max(nlab1);
%     n2 = max(nlab2);
%     n = max(n1,n2); 

	[nlab1,nlab2,lablist] = renumlab(lab1,lab2);
    n = size(lablist,1);
    n1 = size(lablist1, 1);
    n2 = size(lablist2, 1);
    
	% Construct matrix of co-occurences (confusion matrix).

    idx_count = 1;
    error_idx = [];
    error_type_idx = [];
    error_type_lablist = {};
    length_lablist1 = size(lablist1,2);
    length_lablist2 = size(lablist2,2);
    
    lablist1 = char([ cellstr('No_label'); cellstr(lablist)]);
    lablist2 = char([ cellstr('Reject'); cellstr(lablist)]);
    
    
	C = zeros(n+1,n+1);
	for i = 0:n
		K = find(nlab1==i);
		if (isempty(K))
            %Clase no presente
			C(i+1,:) = zeros(1,n+1);
		else
			for j = 0:n
                aux_idx = find(nlab2(K)==j);
				laux_idx = length(aux_idx);
                C(i+1,j+1) = laux_idx;
%                 if(i>0 && j>0 && i~=j && laux_idx > 0)
                if(i~=j && laux_idx > 0)
                    %es un error de un ejemplo etiquetado correctamente.
                    error_idx = [ error_idx; K(aux_idx(:))];
                    error_type_idx = [error_type_idx; repmat(idx_count, laux_idx,1)];                    
                    error_type_lablist{idx_count} = [lablist1(i+1,1:min(length_lablist1,4)) ' x ' lablist2(j+1,1:min(length_lablist2,4))];
                    idx_count = idx_count + 1;
                end
			end
		end
    end
    
    error_type_lablist = char(error_type_lablist);

return
