%MEANCOV Fixed mapping estimating means and covariances of multiclass data
% 
%   [U,G] = MEANCOV(A,N)
%   [U,G] = A*MEANCOV([],N)
%   [U,G] = A*MEANCOV(N)
% 
%  INPUT
%   A	  Dataset
%   N	  Normalization to be useduse for calculating covariances: by M, the 
%	      number of samples in A (N = 1) or by M-1 (default, unbiased, N = 0).
%
% OUTPUT
%   U     Mean vectors, given as a dataset with all annotation of A
%   G     Covariance matrices
%
% DESCRIPTION  
% Computation of a set of mean vectors U and a set of covariance matrices G
% of the C classes in the dataset A. The covariance matrices are stored as a
% 3-dimensional matrix G of the size K x K x C, the class mean vectors as a
% labeled dataset U of the size C x K.
%
% The use of soft labels or target labels is supported.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>) 
% DATASETS, NBAYESC, DISTMAHA

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: meancov.m,v 1.7 2010/02/08 15:31:47 duin Exp $

function [U,G] = meancov(varargin)

	argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],0);
  
  if mapping_task(argin,'definition')
    
    U = define_mapping(argin,'fixed');
    
  elseif mapping_task(argin,'training')			% Compute mapping.
  
    [a,n] = deal(argin{:});

    % N determines whether the covariances are normalized by M (N = 1) or by 
    % M-1 (unbiased, N = 0), where M is the number of objects.

    if (n ~= 1) && (n ~= 0)
      error('Second parameter should be either 0 or 1.')
    end

    if (isdouble(a))			% A is a matrix: compute mean and covariances
      U = mean(a);							% 	in the usual way.
      if nargout > 1
        G = prcov(a,n); 
      end

    elseif (isdatafile(a))

      if nargout  < 2
        U = meancov_datafile(a,n);
      else
        [U,G] = meancov_datafile(a,n);
      end

    elseif (isdataset(a))

      [m,k,c] = getsize(a);
      if nargout == 2
        G = zeros(k,k,c);
      end
      if (islabtype(a,'crisp'))

        if (c==0)  % special solution if all data is unlabeled
          U = mean(+a);
          if nargout > 1
            G = prcov(+a,n); 
          end
        else
          Jall = findclasses(a);
          for i = 1:c     
            %J = findnlab(a,i);
            J = Jall{i};
            if isempty(J)
              U(i,:) = NaN(1,k);
              if nargout > 1
                G(:,:,i) = NaN(k,k);
              end
            else
              U(i,:) = mean(a(J,:),1);
              if (nargout > 1)
                G(:,:,i) = covm(a(J,:),n);	
              end
            end
          end
        end
        labu = getlablist(a);
      elseif (islabtype(a,'soft'))
        problab = gettargets(a);
        % Here we also have to be careful for unlabeled data
        if (c==0)
          prwarning(2,'The dataset has soft labels but no targets defined: using targets 1');
          U = mean(+a);
          if nargout > 1
            G = prcov(+a,n);
          end
        else
          U = zeros(c,k);
          for i = 1:c
            g = problab(:,i);
            if mean(g) == 0 % all targets are zero: zero mean, zero cov.
              G(:,:,i) = zeros(k,k);
            else
              % Calculate relative weights for the means.
              g = g/mean(g); 
              U(i,:) = mean(a.*repmat(g,1,k));	% Weighted mean vectors	
              if (nargout > 1)
                u  = mean(a.*repmat(sqrt(g),1,k));
                % needed to weight cov terms properly
                G(:,:,i) = covm(a.*repmat(sqrt(g),1,k),1) - U(i,:)'*U(i,:) + u'*u;
                % Re-normalise by M-1 if requested.
                if (n == 0)
                  G(:,:,i) = m*G(:,:,i)/(m-1);
                end
              end
            end
          end
        end
        labu = getlablist(a);
      else
        % Default action.
        U = mean(a);
        if nargout > 1
          G = covm(a,n);
        end
        labu = [];
      end

      % Add attributes of A to U.
      U = prdataset(U,labu,'featlab',getfeatlab(a), ...
                                  'featsize',getfeatsize(a));
      if (~islabtype(a,'targets'))
      %	p = getprior(a);
        U = setprior(U,a.prior); 
      end

    else
      error('Illegal datatype')

    end

  end
return



%MEANCOV Datafile overload

function [u,g] = meancov_datafile(a,n)

	if nargin < 2, n = []; end
	
  [m,k,c] = getsize(a);
  k = 0; % just for detecting first loop
	next = 1;
	while next > 0
		[b,next] = readdatafile(a,next);
    if k == 0
      k = size(b,2);
      if c == 0
        u = zeros(1,k);
      else
        u = zeros(c,k);
      end
      if nargout > 1
        if c == 0
          g = zeros(k,k);
        else
          g = zeros(k,k,c);
        end
      end
    end
		%compute for each call to readdatafile contributions to mean and cov
		bb = + b;
		if c == 0
			u = u + sum(bb,1);
		else
			nlab = getnlab(b);
			for j=1:size(b,1)
				if nlab(j) > 0
					u(nlab(j),:) = u(nlab(j),:) + bb(j,:);
					if nargout > 1
						g(:,:,nlab(j)) = g(:,:,nlab(j)) + bb(j,:)'*bb(j,:);
					end
				end
			end
		end
	end

	f = classsizes(a);
	u = u ./ repmat(f',1,k);

	if nargout == 2
		for j=1:c
			g(:,:,j) = g(:,:,j)/f(j) - u(j,:)'*u(j,:);
			if isempty(n) || n == 0
				g(:,:,j) = (f(j)/(f(j)-1))*g(:,:,j);
			end
		end
	end
	u = prdataset(u,getlablist(b));
	u.prior = b.prior;
return
