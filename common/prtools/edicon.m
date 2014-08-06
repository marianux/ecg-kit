% EDICON Multi-edit and condense a training set 
% 
% 	J = EDICON(D,NSETS,NITERS,NTRIES)
%
% INPUT
%   D       Distance matrix dataset
%   NSETS   Number of subsets for editing, or [] for no editing (default: 3)
%   NITERS  Number of iterations for editing (default: 5)
%   NTRIES  Number of tries for condensing, or [] for no condensing (dflt: 10)
%
% OUTPUT
%   J       Indices of retained samples
%
% DESCRIPTION
% Returns the set of objects J such that the nearest neigbour gives zero error 
% on the remaining objects. If MODE = 0, multi-edit the dataset represented by 
% distance matrix D first. D can be computed from a dataset A by A*proxm(A).
% 
% REFERENCES
% Devijver, P. and Kittler, J. "Pattern recognition", Prentice-Hall, 1982.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PRDATASET, KNNC, PROXM

% Copyright: E. Pekalska and D. de Ridder, {E.Pekalska,D.deRidder}@ewi.tudelft.nl
% Faculty of Electrical Engineering, Mathematics and Computer Science 
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function J = edicon (D,no_subsets,max_iter,no_tries)

    if (nargin < 4 | isempty(no_tries))
        prwarning(4,'number of tries for condensing not specified, assuming 10');
        no_tries = 10;
    end;
    if (nargin < 3 | isempty(max_iter))
        prwarning(4,'number of iterations for editing not specified, assuming 5');
        max_iter = 5;
    end;
    if (nargin < 2 | isempty(no_subsets))
        prwarning(4,'number of subsets for editing not specified, assuming 3');
        no_subsets = 5;
    end;

    % Extract dataset information.
    
	nlab    = getnlab(D); lablist = getlablist(D);
	[m,k,c] = getsize(D); p       = getprior(D);

    if (~isdataset(D) | (m ~= k))
        error('require a square distance matrix dataset'); 
    end
    
    J = 1:m;            % Initialise index array.
    
    % If requested, apply the multi-edit algorithm.
    
    if (~isempty(no_subsets))
        iter = 1; 
        while (iter <= max_iter)

      	    % 1. Create NO_SUBSETS distinct subsets out of the sample indices
            %    remaining in J. Store the subset indices in subset_ind{.}.
            %    Note: this is slightly more complicated than strictly
            %    necessary, to enforce that each class is represented equally
            %    in each subset.
        
    		m = length(J); J = J(randperm(m)); 
		    subset_size = floor(m/no_subsets);

            % Find the MC(J) indices of samples belonging to class J in CLASS_IND{J}.
            % Also create a random permutation for class J.
        
            for j = 1:c
			    class_ind{j} = find(nlab(J)==j);
                mc(j)        = length(class_ind{j});
			    perm{j}      = randperm(mc(j));
    		end;
        
            % Distribute the samples over the subsets, per class.
        
    		for i = 1:no_subsets
		        subset_ind{i} = [];
			    for j = 1:c
                    num           = min(floor(mc(j)/no_subsets),length(class_ind{j}));
        		    subset_ind{i} = [subset_ind{i}; class_ind{j}(1:num)];
				    class_ind{j}  = class_ind{j}(num+1:end);
    			end;
			    perm_subset_ind{i} = J(subset_ind{i});
    		end;
    
          	% 2. In turn, classify each subset I using 1-NN on subset I+1.

          	drop = [];
  	        for i = 1:no_subsets
    			L1 = perm_subset_ind{i}; 
                L2 = perm_subset_ind{mod(i,no_subsets)+1};
    			if (length(L2)>1)
	  	    	    [dummy,nearest] = min(D(L2,L1));
    			else
				    nearest = 1;
    			end;
                drop = [drop; subset_ind{i}(find(nlab(L1)~=nlab(L2(nearest))))];
            end;

         	% 3. Discard incorrectly classified samples from J.

  	        if (isempty(drop))
  		        iter = iter + 1;
         	else
   		        J(drop) = []; iter = 0;
          	end;
        end;

        D = D(J,J);
    end;

    % Condense.

    if (~isempty(no_tries))
        % Extract dataset information.
    
    	nlab    = getnlab(D); lablist = getlablist(D);
	    [m,k,c] = getsize(D); p       = getprior(D);

        D = +D;

        % The whole procedure starts from a random object and continues to add 
        % objects, so it is not optimal. Therefore it is repeated NO_TRIES times 
        % and the smallest set found is returned.

        K = zeros(m,no_tries); storesz = zeros(1,no_tries);
        for o = 1:no_tries
        
            % Start with 1 sample index in STORE, the others in GRABBAG.
            p = randperm(m); store = p(1); grabbag = p(2:m);

            % While there are changes...
            transfer = 0;
            while (~transfer) & (~isempty(grabbag)) 
    	        storelab    = nlab(store);
	            grabbaglab  = nlab(grabbag);
	            new_grabbag = []; 

                % For all samples in GRABBAG...
            	transfer = 0;
	            for k = 1:length(grabbag)
                    % ... find the nearest sample in STORE...
   	                [dummy,z] = min(D(grabbag(k),store));
						
                    % ... if it has a different label, move it to STORE.
   	                if (storelab(z) ~= grabbaglab(k))
                        store    = [store; grabbag(k)];
                        storelab = [storelab; grabbaglab(k)];
                        transfer = 1;
               	    else   
                        new_grabbag = [new_grabbag; grabbag(k)];   
                   	end;  
    	        end; 

            	grabbag = new_grabbag;
     
              end

              % Remember the STORE and its size for this repetition.
              storesz(o)        = length(store);
              K(1:storesz(o),o) = store;
  
        end

        % Take the STORE with minimal size.
        [minsz,minszind] = min(storesz);
        J = J(K(1:minsz,minszind));
    end;
    
return
