%FEATEVAL Evaluation of feature set for classification
% 
% 	J = feateval_new(A,CRIT,T)
% 	J = feateval_new(A,CRIT,N)
% 
% INPUT
%       A      input dataset
%       CRIT   string name of a method or untrained mapping
%       T      validation dataset (optional)
%       N      number of cross-validations (optional)
%
% OUTPUT
%       J      scalar criterion value
%
% DESCRIPTION
%  Evaluation of features by the criterion CRIT for classification,
%  using objects in the dataset A. The larger J, the better. Resulting
%  J-values are incomparable over the various methods.
%  The following methods are supported:
%  
%   crit='in-in' : inter-intra distance.
%   crit='maha-s': sum of estimated Mahalanobis distances.
%   crit='maha-m': minimum of estimated Mahalanobis distances.
%   crit='eucl-s': sum of squared Euclidean distances.
%   crit='eucl-m': minimum of squared Euclidean distances.
%   crit='NN'    : 1-Nearest Neighbour leave-one-out
%                  classification performance (default).
%                  (performance = 1 - error). 
%  
%  CRIT can also be any untrained classifier, e.g. LDC([],1e-6,1e-6). 
%  The classification error is used for a performance estimate. If 
%  supplied, the dataset T is used for obtaining an unbiased estimate 
%  of the performance of classifiers trained with the dataset A. 
%  If a number of cross-validations N is supplied, the routine is
%  run for N times with different training and test sets generated
%  from A by cross-validation. Results are averaged. If T nor N are 
%  given, the apparent performance on A is used. 
% 
% SEE ALSO
% DATASETS, FEATSELO, FEATSELB, FEATSELF, FEATSELP, FEATSELM, FEATRANK

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% REVISIONS
% DXD1: David Tax, 08-05-2003
%       I added the inter/intra distance criterion.

% $Id: feateval.m,v 1.8 2005/12/15 20:43:25 duin Exp $

function J = feateval_new(a,crit,t, varargin)

	prtrace(mfilename);
	
	[ma,k,c] = getsize(a);

    
	if nargin < 2
		crit = 'NN';
	end
	if nargin < 3
		t =[]; 
		prwarning(4,'Where needed, input dataset is used for validation')
	end

	if(isscalar(t) && ~isdataset(t) ) % cross-validation desired, t rotations
		
        if( t > 0 )
            K = crossval(a,nmc,t,0); % trick to get rotation set from crossval
            J = 0;
            JALL = 1:size(a,1);
            for j=1:t
                JIN = JALL;
                JOUT = find(K==j);
                JIN(JOUT) = [];
                JOUT = JALL(JOUT);
                train = a(JIN,:);
                test  = a(JOUT,:);
                J = J + feval(mfilename,train,crit,test);
            end
            J = J/t;
            return
        else

            J = cross_val_leave1out_rec_index(a, crit);
            
            return;
            
        end
    end

%     if(~isempty(t))
%         ClassSizes = classsizes(t);
%         bClassesPresent_t = ClassSizes > 0 ;
%     else
%         ClassSizes = classsizes(a);
%         bClassesPresent_t = ClassSizes > 0 ;
%     end
%     
%     if( sum(bClassesPresent_t) < 2 )
%         error('Debe haber representadas al menos 2 clases en el dataset.');
%     end
    
% 	islabtype(a,'crisp');
	
    isvaldset(a,1,2); % at least 1 object per class, 2 classes
	iscomdset(a,t);


    
	if ischar(crit)
		%DXD1
		if strcmp(crit,'in-in')     % inter/intra distances
			islabtype(a,'crisp','soft');
            if isempty(t)
				[U,G] = meancov_new(a,0);
                prior = getprior(a);
			else
 				[U,G] = meancov_new(t,0);
                prior = 1/sum(classsizes(t))*classsizes(t);
%                 S_w = sum(G(:,:,bClassesPresent_t),3);
%                 
%                 iClassMean = +U;
%                 iTotalMean = prior(bClassesPresent_t) * iClassMean(bClassesPresent_t,:);
%                 
%                 S_b = zeros( size(iClassMean,2) );
%                 
%                 for j = find(bClassesPresent_t)
%                     iAux = iClassMean(j,:) - iTotalMean; % between scatter
%                     S_b = S_b + prior(j)*(iAux' * iAux);
%                 end
                
            end

            S_w = reshape(nansum(reshape(G(:,:,bClassesPresent_t),k*k,sum(bClassesPresent_t))*prior(bClassesPresent_t)',2),k,k); % within scatter
            [dummy, S_b] = meancov_new(+U, 0); % between scatter
            
            if( det(S_w) == 0 )
                
                J = 0;
                
            else
                
                J = trace(inv(S_w)*S_b);
            
            end

%             t
%             U
%             +U
%             disp(['G:' num2str(G) ' Sw:' num2str(S_w) ' S_b:' num2str(S_b)]);
            
        elseif( strcmp(crit,'maha-s') || strcmp(crit,'maha-m') )% Mahalanobis distances
			islabtype(a,'crisp','soft');
			if isempty(t)
				D = distmaha_new(a,[],[]);
			else
				[U,G] = meancov_new( a, 0 );
				D = distmaha(t,U,G);
				D = meancov_new(D,0);
            end
            D = +D;
            D = D(bClassesPresent_t,bClassesPresent_t);

            if strcmp(crit,'maha-m')
				D = D + realmax*eye(sum(bClassesPresent_t));
				J = min(min(D));
			else
				J = sum(sum(D))/2; 
            end
            
		elseif strcmp(crit,'eucl-s') || strcmp(crit,'eucl-m') % Euclidean distances
			islabtype(a,'crisp','soft');
			U = meancov_new(a,0);
			if(isempty(t))
				D = distm_new(U,U);
            else
				D = distm_new(t,U);
				D = meancov_new(D,0);
            end
            D = +D;
            D = D(bClassesPresent_t,bClassesPresent_t);

			if strcmp(crit,'eucl-m')
				D = D + realmax*eye( sum(bClassesPresent_t) );
				J = min(min(D));
			else
				J = sum(sum(D))/2; 
			end
		elseif strcmp(crit,'NN')	% 1-NN performance
			islabtype(a,'crisp','soft');
			if isempty(t)
				J = 1 - testk(a,1);
			else
				J = 1 - testk(a,1,t);
			end
		else
			error('Criterion undefined');
		end
	else
		ismapping(crit);
		isuntrained(crit);
		if isempty(t)
			J = 1 - (a * (a * crit) * testc);
		else
			J = 1 - (t * (a * crit) * testc);
		end
	end

return
