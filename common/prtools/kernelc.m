%KERNELC Arbitrary kernel/dissimilarity based classifier
%
%   W = KERNELC(A,KERNEL,CLASSF)
%   W = A*KERNELC([],KERNEL,CLASSF)
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
% A = GENDATB([100 100]);    % Training set of 200 objects
% R = GENDATB([10 10]);      % Representation set of 20 objects
% V = KERNELM(R,'p',3);      % Compute kernel
% W = KERNELC(A,V,FISHERC)   % compute classifier
% SCATTERD(A);               % Scatterplot of trainingset
% HOLD ON; SCATTERD(R,'ko'); % Add representation set to scatterplot
% PLOTC(W);                  % Plot classifier
%
% SEE ALSO
% DATASETS, MAPPINGS, KERNELM, PROXM, USERKERNEL

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = kernelc(a,kernel,classf)

		
	if nargin < 3 | isempty(classf), classf = fisherc; end
	if nargin < 2 | isempty(kernel)
		kernel = kernelm([],[],'random',0.7,100); 
	end
	
	if nargin < 1 | isempty(a)
		w = prmapping(mfilename,'untrained',{kernel,classf}); 
	else % training
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

	

