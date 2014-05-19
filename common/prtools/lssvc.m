function W = lssvc(A, TYPE, PAR, C)
%LSSVC Least-Squares Support Vector Classifier
%
%       W = lssvc(A,TYPE,PAR,C);
%
% INPUT
%   A       dataset
%   TYPE    Type of the kernel (optional; default: 'p')
%   PAR     Kernel parameter (optional; default: 1)
%   C       Regularization parameter (optional; default: 1)
%
% OUTPUT
%   W       Mapping: Least-Squares Support Vector Classifier
% 
% DESCRIPTION
% Optimizes a least-squares support vector classifier for the dataset A by 
% quadratic programming. The classifier can be of one of the types 
% as defined by PROXM. Default is linear (TYPE = 'p', PAR = 1). The 
% regularization parameter C allows for seeting the level of overfitting. 
% A smaller value for C allows for more class overlap. Default C = 10.
%
% NOTE
% This implementation uses the LS-SVMLab toolbox. This toolbox is available
% from http://www.esat.kuleuven.ac.be/sista/lssvmlab/. Make sure the path
% to LS-SVMLab is set in Matlab.
% 
% See also MAPPINGS, DATASETS, PROXM

% Copyright: L.J.P. van der Maaten, l.vandermaaten@micc.unimaas.nl
% MICC-IKAT, Maastricht University, Maastricht, The Netherlands

    name = 'LSSVC';
    
    % Check whether LS-SVMLab is installed and set in the path
    if ~exist('trainlssvm.m', 'file')
        error('Could not find LS-SVMLab. Make sure LS-SVMLab is installed and set in the Matlab path.');
    end
 
    % Perform LS-SVM training on dataset A
    if nargin <= 1 || (nargin >= 2 && ~isa(TYPE, 'prmapping'))
        
        % Set LS-SVM parameters
        if ~exist('TYPE', 'var') || ~isstruct(TYPE)
            if nargin < 2 || isempty(TYPE), data.kernel = 'p';    else data.kernel = TYPE; end       % kernel function
            if nargin < 3 || isempty(PAR),  data.par = 1;         else data.par = PAR; end           % kernel parameters           
            if nargin < 4 || isempty(C),    data.gamma = 10;      else data.gamma = C; end           % regularization parameter
            data.par2 = -1;

            % LS-SVMLab uses different naming for kernels; correct this
            switch lower(data.kernel)
                case {'p', 'polynomial'},       data.kernel = 'poly_kernel'; data.par2 = 1;
                case {'h', 'homogeneous'},      data.kernel = 'poly_kernel'; data.par2 = 0;
                case {'e', 'exponential'},      error('This kernel is not supported by LS-SVMLab.');
                case {'r', 'radial_basis'},     data.kernel = 'RBF_kernel'; data.par = data.par^2;          % LS-SVMLab uses sigma^2 as input
                case {'s', 'sigmoid'},          data.kernel = 'MLP_kernel';
                case {'d', 'distance'},         error('This kernel is not supported by LS-SVMLab.');
                case {'m', 'minkowski'},        error('This kernel is not supported by LS-SVMLab.');
                case {'c', 'city-block'},       error('This kernel is not supported by LS-SVMLab.');
                case {'o', 'cosine'},           error('This kernel is not supported by LS-SVMLab.');
                otherwise, error('Unknown kernel function.');
            end
        else
            data = TYPE;
        end
        
        % Handle the case in which the dataset is empty
        if ~exist('A', 'var') || isempty(A)
            W = prmapping('lssvc', 'untrained', data); 
            W = setname(W, name);
            return;
        end
        
        % Perform some checks on the data
        islabtype(A, 'crisp');          % allow crisp labels only 
        isvaldset(A, 1, 2);             % at least one object per class, two objects
        if ischar(A.lablist), error('Only numerical labels are allowed in LSSVC.'); end
        
        % Train binary LS-SVM
        [m, k, c] = getsize(A);         % size of the training set
		if length(A.lablist) <= 2
            if strcmp(data.kernel, 'poly_kernel')
                [alpha, b] = trainlssvm({A.data, A.labels, 'classification', data.gamma, [data.par2; data.par], data.kernel});
            else
                [alpha, b] = trainlssvm({A.data, A.labels, 'classification', data.gamma, data.par, data.kernel});
            end
		
		% Train multi-class LS-SVM
		else
			[labels_code, codebook, old_codebook] = code(A.labels, 'code_MOC');
            if strcmp(data.kernel, 'poly_kernel')
                [alpha, b] = trainlssvm({A.data, labels_code, 'classification', data.gamma, [data.par2; data.par], data.kernel});
            else
                [alpha, b] = trainlssvm({A.data, labels_code, 'classification', data.gamma, data.par, data.kernel});
            end
		end

        % Store results of training in data-struct
        data.A = A;
        data.alpha = alpha;
        data.b = b;
		if length(A.lablist) > 2
			data.labels_code = labels_code;
			data.codebook = codebook;
			data.old_codebook = old_codebook;
		end
        W = prmapping('lssvc', 'trained', data, getlablist(A), k, c); 
        W = setname(W, name);
    
        
    % Apply LS-SVM mapping B on dataset A
    elseif nargin == 2 && isa(TYPE, 'prmapping')
    
        % Initialize some variables
        B = TYPE;
        [m, k] = getsize(A);        % size of the test set
        [k, c] = size(B);           % K features with C classes
        W = zeros(m, c);            % output: C class densities for M objects
        
        % Perform binary classification 
		if length(B.data.A.lablist) <= 2
            if strcmp(B.data.kernel, 'poly_kernel')
                C = simlssvm({B.data.A.data, B.data.A.labels, 'classification', B.data.gamma, [B.data.par2; B.data.par], B.data.kernel, 'preprocess'}, ...
                             {B.data.alpha, B.data.b}, A.data);
            else
				C = simlssvm({B.data.A.data, B.data.A.labels, 'classification', B.data.gamma, B.data.par, B.data.kernel, 'preprocess'}, ...
                             {B.data.alpha, B.data.b}, A.data);
            end
            
		% Perform multi-class classification (only numerical labels!)
        else
            if strcmp(B.data.kernel, 'poly_kernel')
                C = simlssvm({B.data.A.data, B.data.labels_code, 'classification', B.data.gamma, [B.data.par2; B.data.par], B.data.kernel, 'preprocess'}, ...
                             {B.data.alpha, B.data.b}, A.data);
            else
                C = simlssvm({B.data.A.data, B.data.labels_code, 'classification', B.data.gamma, B.data.par, B.data.kernel, 'preprocess'}, ...
                             {B.data.alpha, B.data.b}, A.data);
            end
			C = code(C, B.data.old_codebook, [], B.data.codebook, 'codedist_hamming');
        end
                                
        % Make sure all labels are numbered in the same way as in A
        for i=1:length(C)
            [tf, loc] = ismember(B.data.A.lablist, C(i));
            W(i, find(loc)) = 1;
        end
        
        % Return new labeled dataset
        W = setdat(A, W, B);
        
        
    % Should not happen
    else
        error('Illegal call.');
    end
