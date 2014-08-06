%ROC Receiver-Operator Curve, deprecated
% 
%   E = ROC(A,W,C,N)
%   E = ROC(B,C,N)
%
% INPUT
%   A  Dataset
%   W  Trained classifier, or
%   B  Classification result, B = A*W*CLASSC
%   C  Index of desired class (default: C = 1)
%   N  Number of points on the Receiver-Operator Curve (default: 100)
%
% OUTPUT
%   E  Structure containing the error of the two classes
%
% DESCRIPTION
% Computes N points on the Receiver-Operator Curve (ROC)of the classifier W
% for class C in the labeled dataset B, which is typically the result of
% B = A*W; or for the dataset A labelled by applying the (cell array of)
% trained classifier(s) W.
%
% Note that a ROC is related to a specific class (class C) for which the
% errors are plotted horizontally. The total error on all other classes is
% plotted vertically. The class index C refers to its position in the label
% list of the dataset (A or B). It can be found by GETCLASSI.
%
% The curve is computed for N thresholds of the posteriori probabilities
% stored in B. The resulting error frequencies for the two classes are
% stored in the structure E. E.XVALUES contains the errors in the first
% class, E.ERROR contains the errors in the second class. In multi-class
% problems these are the mean values in a single class, respectively the
% mean values in all other classes. This may not be very useful, but not
% much more can be done as for multi-class cases the ROC is equivalent to a
% multi-dimensional surface.
%
% Use PLOTE(E) for plotting the result. In the plot the two types of error
% are annotated as 'Error I' (error of the first kind) and 'Error II' (error
% of the second kind). All error estimates are weighted according the class
% prior probabilities. Remove the priors in A or B (by setprior(A,[])) to
% produce a vanilla ROC.
%
% This routine calls PRROC and returns its result.
%
% EXAMPLES
%	Train set A and test set T:
%	  B = T*NMC(A); E = PRROC(T,50); PLOTE(E); % Plots a single curve
%	  E = PRROC(T,A*{NMC,UDC,QDC});  PLOTE(E); % Plots 3 curves
%
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd edition, 
%    John Wiley and Sons, New York, 2001.
% 2. A. Webb, Statistical Pattern Recognition, John Wiley & Sons, New York, 
%    2002.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PLOTE, REJECT, TESTC, GETCLASSI, PRROC
