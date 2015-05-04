%EMCLUST Expectation-Maximization clustering
%
%  [LABELS,W_EM] = EMCLUST (A,W_CLUST,K,LABTYPE,FID)
%
% INPUT
%   A         Dataset, possibly labeled
%   W_CLUST   Cluster model mapping, untrained (default: nmc)
%   K         Number of clusters (default: 2)
%   LABTYPE   Label type: 'crisp' or 'soft' (default: label type of A)
%   FID       File ID to write progress to (default [], see PRPROGRESS)
%
% OUTPUT
%   LABELS    Integer labels for the objects in A pointing to their cluster
%   W_EM      EM clustering mapping
%
% DESCRIPTION
% The untrained classifier mapping W_CLUST is used to update an initially
% labeled dataset A by iterating the following two steps:
%   1. Train W   :  W_EM = A*W_CLUST
%   2. Relabel A :  A    = dataset(A,labeld(A*W_EM*classc))
% This is repeated until the labeling does not change anymore. The final
% classification matrix is returned in B. The final crisp labeling is returned
% in LABELS. W_EM may be used for assigning new objects.
%
% If K is given, a random initialisation for K clusters is made and labels
% of A are neglected. 
%
% LABTYPE determines the type of labeling: 'crisp' or 'soft'. Default: label
% type of A. It is assumed W_CLUST can handle the LABTYPE requested.
% Only in case LABTYPE is 'soft' the traditional EM algorithm is followed.
% In case LABTYPE is 'crisp' EMCLUST follows a generalised k-means
% algorithm.
%
% SEE ALSO
% MAPPINGS, DATASETS, KMEANS, PRPROGRESS

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: emclust.m,v 1.9 2009/02/03 21:07:26 duin Exp $

function [new_lab,w_em] = emclust (a,w_clust,n,type,fid)

	prtrace(mfilename);

	n_ini		= 500;			% Maximum size of subset to use for initialisation.
	epsilon = 1e-6;			% Stop when average labeling change drops below this.

	% Check arguments.
	if (nargin < 5), fid = []; end
	if (nargin < 4)
		prwarning(3,'No label type specified, using label type of dataset A.');
		type = []; 
	end
	if (nargin < 3) | isempty(n)
		prwarning(3,'No number of clusters specified, using number of classes in A.');
		n = []; 
	end
	if (nargin < 2) | isempty(w_clust)
		prwarning(2,'No clustering mapping specified, assuming NMC.');
		w_clust = nmc;   
	end

  isuntrained(w_clust);   % Assert that clustering mapping is untrained.

  % Determine number of clusters N and initialisation method.

	a = testdatasize(a); 
	islabtype(a,'crisp','soft');
	[m,k,c] = getsize(a); 
	rand_init = 1;
	if (isempty(n))
		if (c == 1)						% For one class, find two clusters.
			n = 2;
		else
			n = c;							
			rand_init = 0; 			% Use given classification as initialisation.
		end
	end

	if (n < 1),  error('Number of clusters should be at least one.'); end
	if (n == 1), prwarning(4,'Clustering with 1 cluster is trivial.'); end

	% Set label type, if given.

	if ~isempty(type), a = setlabtype(a,type); end
	a = setprior(a,[]); % make sure that priors will be deleted
	
	% Initialise by performing KCENTRES on...

	prwaitbar(2,'EM Clustering, initialization');
	prwaitbar(2,1);
	if (rand_init)
		
        not_found = 1;
		itern = 0;
		while(not_found)
			% try to find an initialisation with all class sizes > 1
            
            if (m > n_ini)       						% ... a random subset of A.
                prwarning(2,'Initializing by performing KCENTRES on a subset of %d samples.', n_ini);
    % 			a_ini = +gendat(+a,n_ini);	
                %sampleo uniformemente el feature-space
                a_ini = +a;
                a_ini = sum(abs(a_ini),2);
                [~, sort_idx ] = sort(a_ini);
                a_ini = +a(sort_idx(randsample(m, n_ini)),:);
            else
                prwarning(2,'Initializing by performing KCENTRES on the training set.');
                a_ini = +a;								% ... the entire set A.
            end
            
			itern = itern + 1;
			if itern > 3
				error('Not possible to find desired number of components')
			end
			% add some noise to data to avoid problems
			% 50 trials
            assign = zeros(1,n-1);
            std_a_ini = std(a_ini);
            jj = 1;
            iter_max = 50;
            
            while( jj <= iter_max )
                try
                    assign  = kcentres(+distm(a_ini.*(ones(size(a_ini))+ bsxfun(@times, max( [repmat(0.001,size(std_a_ini)); 0.01*std_a_ini] ), randn(size(a_ini))))),n,50);
                    if( length(unique(assign)) == n )
                        a_ini = prdataset(a_ini,assign); 
                        a_ini = setprior(a_ini,getprior(a_ini,0));
                        if( isvaldfile(a_ini,1,2) )
                            break 
                        end
                    end

                catch ME
%                     if( strcmpi(ME.message, 'kcentres fails as some objects are identical: add some noise'))
%                         std_a_ini = std_a_ini * 2;
%                     else
%                         rethrow(ME)
%                     end
                    error('Not possible to find desired number of components')
                end
                
                jj = jj + 1;
                
            end
            
            if( jj > iter_max )
                error('Not possible to find desired number of components')
            end

            
            % Train initial classifier on labels generated by KCENTRES and find
            % initial hard labels. Use NMC instead of W_CLUST to make sure that we 
            % always have enough data to estimate the parameters.
            
            try
                d = a*(a_ini*nmc);
            catch ME
                fprintf(2, 'Valid %d\n\n', isvaldfile(a_ini,1,2) )  
                error('Not possible to find desired number of components')
            end
            if (islabtype(a,'soft'))
				new_lab = +d;
				not_found = 0;
			else
				new_lab = d*labeld;
                cs_new_lab = classsizes(prdataset(d,new_lab));
				if length(cs_new_lab) == n && all( cs_new_lab > 1)
					not_found = 0;
				end
			end
		end
		lablist_org = [];
	else
		lablist_org = getlablist(a);
		a = setlablist(a,[1:c]');
		new_lab = getlabels(a);		% Use given labeling.
	end

	% Ready for the work.
	iter = 0;
	change = 1;
	prwaitbar(2,2,'EM Clustering, EM loop')
	prwaitbar(100,['using ' getname(w_clust)]);
  if (islabtype(a,'soft'))
		prprogress(fid,'emclust optim: iter, change (mse):\n');
		prprogress(fid,' %i, %f \n',0,0);
		a = setlabels(a,new_lab);
		a = setprior(a,getprior(a,0));
		laba = getlabels(a);
		lab = new_lab;
  	while (change > epsilon)       	% EM loop, run until labeling is stable.
			prwaitbar(100,100-100*exp(-iter/10));
  		w_em = a*w_clust;             % 1. Train classifier, density output.
  		b = a*(w_em*classc);          % 2. Assign probability to training samples.
  		a = settargets(a,b);          % 3. Insert probabilities as new labels.
  		change = mean(mean((+b-lab).^2)); lab = b;          
			prprogress(fid,' %i, %f \n',iter,change);
			iter = iter+1;
			if iter > 500
				prwarning(1,'emclust stopped after 500 iterations')
				change = 0;
			end
		end
	else  % crisp labels
		
        prprogress(fid,'emclust optim: iter, change (#obj), #clust:\n');
		prprogress(fid,' %i, %i %i \n',0,0,0);
		lab = ones(m,1);
        
        while (change ~= 0 && any(lab ~= new_lab)) % EM loop, run until labeling is stable.
            
            prwaitbar(100,100-100*exp(-iter/10));
            a = setlabels(a,new_lab); 		% 0. Set labels and store old labels.
            a = setprior(a,getprior(a,0));%    Set priors to class frequencies
            lab = new_lab;								% 
            a = remclass(a,1);            %    demand class sizes > 2 objects
            itern = 0;
            
            while getsize(a,3) < n        %    increase number of classes if necessary
                
                itern = itern + 1;
                
                if itern > 5
                    error('Not possible to find desired number of components')
                end
                laba = getlablist(a);
                labmax = max(laba);
                N = classsizes(a);
                [Nmax,cmax] = max(N);        % find largest class
                aa = seldat(a,cmax);         % select just that one

                std_aa = std(+aa);
                bContinue = true;
                while( bContinue )
                    try
                        new_lab_aa = kmeans(aa .* (ones(size(+aa))+ bsxfun(@times, max( [repmat(0.001,size(+aa)); 0.01*std_aa] ), randn(size(+aa)))) ,2);   % split it by kmeans
                        bContinue = false;
                    catch ME
                        %                         if( strcmpi(ME.message, 'kcentres fails as some objects are identical: add some noise'))
                        %                             std_aa = std_aa * 2;
                        %                         else
                        %                             rethrow(ME)
                        %                         end
                        error('Not possible to find desired number of components')
                    end

                end



                N1 = sum(new_lab_aa == 1);   
                N2 = sum(new_lab_aa == 2);
                if (N1 > 1 & N2 > 1) % use it if both classes have more than one sample
                    J = findlabels(a,laba(cmax,:));
                    a = setlabels(a,new_lab_aa + labmax,J);
                end
            end

            std_a = std(+a);
            w_em = (a .* (ones(size(+a))+ bsxfun(@times, max( [repmat(0.001,1,size(+a,2)); 1e-6*std_a] ), randn(size(+a))))) * w_clust;             % 1. Compute classifier, crisp output.
            b = a*w_em;               		% 2. Classify training samples.
            new_lab = labeld(b);      		% 3. Insert classification as new labels.
            prprogress(fid,' %i, %i %i \n', ...
            iter,length(find(lab ~= new_lab)),length(unique(new_lab)));
            iter = iter+1;             %DXD Added also the iter for the crisp labels
            if iter > 50
                prwarning(1,'emclust stopped after 50 iterations')
                change = 0;
            end
            
            %para ver la velocidad de convergencia.
%             disp(sum(lab ~= new_lab))
            
        end
    end
	prwaitbar(0)
	prwaitbar(0)
	
% 	if ~isempty(lablist_org) % substitute original labels if desired
% 		new_lab = lablist_org(new_lab);
%         confmat_new
% 		wlab = getlabels(w_em);
% 		wlab = lablist_org(wlab);
% 		w_em = setlabels(w_em,wlab);
% 	end
		
return;
