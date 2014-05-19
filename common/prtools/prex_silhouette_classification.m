%PREX_SILHOUETTE_CLASSIFICATION PRTools introductory example
%
% Presents the construction of a dataset from a set of images,
% builds a classifier and performs an evaluation
%
help prex_silhouette_classification

delfigs
echo on
	% load the Kimia image collection
	a = kimia;
	show(a,18);
	%
	% Select 5 classes and show them
	b = seldat(a,[1 3 5 7,13]);
	figure; show(b);
	%
	% compute the area
	area = im_stat(b,'sum');
	%
	% compute the perimeter
	bx = abs(filtim(b,'conv2',{[-1 1],'same'}));
	by = abs(filtim(b,'conv2',{[-1 1]','same'}));
	bor = or(bx,by);
	figure; show(bor);
	perimeter = im_stat(bor,'sum');
	%
	% construct a dataset with equal class priors and label the features
	c = prdataset([area perimeter]);
	c = setprior(c,0);
	c = setfeatlab(c,char('area','perimeter'));
	%
	% Show the 2D scatterplot
	figure; scatterd(c,'legend');
	%
	% Compute a linear classifier and show it in the scatterplot
	w = ldc(c);
	plotc(w,'col');
	showfigs
	%
	% Compute and print classification errors
	train_error = c*w*testc;
	test_error = crossval(c,ldc,2);
	
echo off
	
	fprintf('\nerror in training set:        %5.3f\n',train_error)
	fprintf('2-fold crossvalidation error: %5.3f\n',test_error)
		
