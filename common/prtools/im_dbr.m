%IM_DBR Image Database Retrieval GUI
%
%    [RANK,TARG,OUTL] = IM_DBR(DBASE,FSETS,CLASSF,COMB)
%
% INPUT
%   DBASE   - Dataset or datafile with N object images
%   FSETS   - Cell array with maximum 4 feature sets
%   CLASSF  - Cell array with untrained classifiers (Default: KNNC([],1))
%   COMB    - Combining classifier (Default: MEANC)
%
% OUTPUT
%   RANK    - Index array ranking the N object images
%   TARG    - Index array pointing to user defined target images
%   OUTL    - Index array pointing to user defined outlier images
%
% DESCRIPTION
% This command generates a Graphical User Interface (GUI) enabling the user
% to label a database of images in 'target' and 'outlier' images in an
% interactive and iterative way. Up to four feature sets can be given and
% corresponding classifiers that assist the user by predict an object ranking 
% based on classification confidences for the 'target' class.
%
% The GUI shows the top-10 of the ranking and the user should classify
% them as targets or outliers (original object labels in DBASE are
% neglected). There are buttons for browsing through the ranked database
% or through the selected targets and outliers. Classifiers can be trained 
% according to two different strategies using the top right buttons:
% Classify - uses all stored target and outlier objects (shown in the top
%            left windows) for building a training set as well as the
%            hand labeled images in the present screen.
% Label    - uses just the hand labeled images in the present screen
%            and neglects the stored targets and outliers. This enables 
%            a more flexible, but still controlled browsing throug the
%            database.
% Reset    - Resets the entire procedure by deleting all selected targets
%            and outliers.
% Quit     - Deletes the GUI and returns the ranking and selected targets
%            and outliers to the user.
% A few additional buttons and sliders for controlling the system behavior:
% - Delete and move buttons for the selected targets and outliers
% - Weights for the feature sets. For each feature set a different
%   classifier is computed generating target confidences for all images. 
%   This influences the operation of the combiniong classifier.
%   The weights can be changed by a slider for every feature set.
%   By default weights are 1.
% - Two buttons for setting all labels as target ('All target') or outlier
%   ('All outlier').
% - Labels for the individual images can be changed by a mouse-click in the
%   image or on the image check-box.
% - For all images a target confidence is computed. Depending on the 'all'
%   and 'unlabeled' radio buttons at the bottom the ranking of all images
%   or of the yet unlabeled images are shown.
% Note: It is not an error, but for most classifiers useless or
% counterproductive to label an object as target as well as outlier.
%
% EXAMPLE
% % This example assumes that the Kimia images are available as datafile
% % and that the DipImage image processing package is available.
% prwaitbar on
% a = kimia_images;
% x = im_moments(a,'hu');
% x = setname(x,'Hu moments');
% y = im_measure(a,a,{'size','perimeter','ccbendingenergy'});
% y = setname(y,'Shape features');
% [R,T,L] = im_dbr(a,{x,y});  % do your own search
% delfigs
% figure(1); show(a(R,:)); % show ranking
% figure(2); show(a(T,:)); % show targets
% figure(3); show(a(L,:)); % show outliers
% showfigs
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, MAPPINGS, KNNC, MEANC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [R,T,L] = im_dbr(dbase,featsets,classf,comb);

	if sscanf(version('-release'),'%i') < 14
		error('IM_DBR needs Matlab version 14 or higher')
	end

	if nargin < 4, comb = meanc; end
	if nargin < 3, classf = knnc([],1); end
	[R,T,L] = image_dbr(dbase,featsets,classf,comb);

return
