%MATCHLAB Compare two labellings and rotate the labels for an optimal match
%
%  LABELS = MATCHLAB(LAB1,LAB2)
%
% INPUT
%   LAB1,LAB2  Label lists of the same objects   
%
% OUTPUT
%   LABELS     A rotated version of LAB2, optimally matched with LAB1
%
% DESCRIPTION
% LAB1 and LAB2 are label lists for the same objects. The returned LABELS 
% constitute a rotated version of LAB2, such that the difference with 
% LAB1 is minimized. This can be used for the alignment of cluster results.
%
% EXAMPLES
% See PREX_MATCHLAB.
%
% SEE ALSO
% DATASETS, HCLUST, MODESEEK, KMEANS, EMCLUST

% $Id: matchlab.m,v 1.2 2006/03/08 22:06:58 duin Exp $


function lab = matchlab_new(lab1,lab2)

	% Compute the confusion matrix and renumber the labels.
	[C, ~, lablist1, lablist2] = confmat(colvec(lab1), colvec(lab2) ); 			

    [~, relab_idx ] = max(C, [], 1);
    
    lab = nan(length(lab2),1);
    for ii = 1:length(relab_idx)
        lab( lab2 == lablist2(ii) ) = lablist1(relab_idx(ii)); 
    end
    
return;

