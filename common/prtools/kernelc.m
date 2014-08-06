%KERNELC Arbitrary kernel/dissimilarity based classifier
%
%   W = KERNELC(A,KERNEL,CLASSF)
%   W = A*KERNELC([],KERNEL,CLASSF)
%   W = A*KERNELC(KERNEL,CLASSF)
%
% INPUT
%   A       Dateset used for training
%   KERNEL  - untrained mapping to compute kernel by A*(A*KERNEL) for
%             training CLASSF or B*(A*KERNEL) for testing with dataset B, or
%           - trained mapping to compute a kernel A*KERNEL for training
%             CLASSF or B*KERNEL for testing with a dataset B
%           KERNEL should be a functions like PROXM, KERNELM or USERKERNEL.
%           Default: reduced dissimilarity representation:
%           KERNELM([],[],'random',0.1,100);
%   CLASSF  Classifier used in kernel space, default LOGLC.
%
% OUTPUT
%   W       Resulting, trained classifier
%
% DESCRIPTION
% This routine defines a classifier W in the input feature space based 
% on a kernel or dissimilarity representation defined by KERNEL and a
% classifier CLASSF to be trained in the kernel space.
%
% In case KERNEL is defined by V = KERNELM( ... ) this routine
% is identical to W = A*(V*CLASSF), Note that if KERNEL is a mapping, it 
% may be trained as well as untrained. In the latter case A is used to
% build the kernel space as well as to optimize the classifier (like in SVC).
%
% EXAMPLE
% a = gendatb([100 100]);    % training set of 200 objects
% r = gendatb([10 10]);      % representation set of 20 objects
% v = proxm(r,'p',3);        % compute kernel
% w = kernelc(a,v,fisherc)   % compute classifier
% scatterd(a);               % scatterplot of trainingset
% hold on; scatterd(r,'ko'); % add representation set to scatterplot
% plotc(w);                  % plot classifier
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, KERNELM, PROXM, USERKERNEL

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = kernelc(varargin)

	argin = shiftargin(varargin,'prmapping');
  argin = setdefaults(argin,[],kernelm([],[],'random',0.7,100),fisherc);
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained');
  else % training
    [a,kernel,classf] = deal(argin{:}); 
		islabtype(a,'crisp');
		isvaldfile(a,1,2); % at least 1 object per class, 2 classes
		a = testdatasize(a,'objects');
		isuntrained(classf);
		ismapping(kernel);
		if isuntrained(kernel)
			kernel = a*kernel;
		end
		K = a*kernel;
		w = kernel*(K*classf);
	end
	w = setname(w,'Kernel Classifier');
	
return

	

