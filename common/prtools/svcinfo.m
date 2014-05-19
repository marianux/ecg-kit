%SVCINFO More information on Support Vector Classifiers
% 
%   [W,J,C,NU,ALGINF] = SVC(A,KERNEL,C,OPTIONS)
%    W                = A*SVC([],KERNEL,C,OPTIONS)
%   [W,J,NU,C,ALGINF] = NUSVC(A,KERNEL,NU,OPTIONS)
%    W                = A*SVC([],KERNEL,NU,OPTIONS)
%
% INPUT
%   A	      Dataset
%   KERNEL  - Untrained mapping to compute kernel by A*(A*KERNEL) during
%             training, or B*(A*KERNEL) during testing with dataset B.
%           - string to compute kernel matrices by FEVAL(KERNEL,B,A)
%           Default: linear kernel (KERNELM([],'p',1));
%   C       Regularization parameter (optional; default: 1)
%   NU      Regularization parameter (0 < NU < 1): expected fraction of SV 
%           (optional; default: max(leave-one-out 1_NN error,0.01))
%   OPTIONS Additional options, see below.
%
% OUTPUT
%   W       Mapping: Support Vector Classifier
%   J       Object indices of support objects		
%   C       Regularization parameter which gives the same classifier by SVC
%   NU      NU parameter which gives the same classifier by NUSVC
%   ALGINF  Structure with additional training information
%
% DESCRIPTION
% Optimizes a support vector classifier for the dataset A by quadratic
% programming. The non-linearity is determined by the kernel.
% If KERNEL = 0 it is assumed that A is already the kernelmatrix (square).
% In this case also a kernel matrix should be supplied at evaluation by B*W 
% or PRMAP(B,W).
%
% If C or NU is NaN this regularisation parameter is optimised by REGOPTC.
%
% The quadratic optimisation is controlled by routines SVO and NUSVO. They make use 
% of one of the following routines, if available:
% - QLD.DLL (Windows) or QLD.MEXxxx under Linux
% - QUADPROG.M in Matlab's optimisation toolbox
% - Matlab's QP.M
%
% The following options are available for fine-tuning the SVC routines
% OPTIONS
%    .MEAN_CENTRING     subtract data mean before the kernel computation (default: 1)
%    .PD_CHECK          force positive definiteness of the kernel by adding a small constant 
%                       to a kernel diagonal (default: 1)
%    .BIAS_IN_ADMREG    it may happen that bias of svc (b term) is not defined, then 
%                       if BIAS_IN_ADMREG == 1, b will be pu in the midpoint of its admissible 
%                       region, otherwise (BIAS_IN_ADMREG == 0) the situation will be considered 
%                       as an optimization failure and treated accordingly (default: 1)
%    .ALLOW_UB_BIAS_ADMREG  (NUSVC only)
%                           it may happen that bias admissible region is unbounded; 
%                           if ALLOW_UB_BIAS_ADMREG == 1, b will be heuristically taken 
%                           from its admissible region, otherwise (ALLOW_UB_BIAS_ADMREG == 0) 
%                           the situation will be considered as an optimization failure and 
%                           treated accordingly (default: 1)
%    .PF_ON_FAILURE     if optimization failed (or bias is undefined and BIAS_IN_ADMREG is 0)
%                       and PF_ON_FAILURE == 1, then Pseudo Fisher classifier will be computed, 
%                       otherwise (PF_ON_FAILURE == 0) an error will be issued (default: 1)
%    .MULTICLASS_MODE   if the multiclass problem has to be solved, MULTICLASS_MODE defines
%                       how it is going to be split in 2-class subproplems: 'single' means
%                       one-against-the rest and 'multi' means
%                       one-against-one (default: 'single')

