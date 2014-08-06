function W = som(a,k,nrruns,eta,h)
%SOM Simple routine computing a Self-Organizing Map (SOM)
%
%           W =  SOM(X,K)
%
% Train a 2D SOM on dataset X. In K the size of the map is defined. The
% map can maximally be 2D. When K contains just a single value, it is
% assumed that a 1D map should be trained. The output of the mapping
% contains the (negative) distances to all neurons. To obtain the index
% of the closest neuron, do A*W*LABELD.
%
%           W =  SOM(X,K,NRRUNS,ETA,H)
%
% Train the SOM for NRRUNS iterations, using learning rate ETA and a
% Gaussian neighborhood function with width H*sqrt(MAXD), where MAXD
% is the maximum distance in the dataset X.
%
% There is the extra feature, that NRRUNS, ETA and H can be vectors,
% such that it can be run several iterations using larger ETA and H,
% and after that with smaller values.
%
% Default: K=[5 5], NRRUNS = [20 40 40], ETA = [0.5 0.3 0.1],
% H = [0.6 0.2 0.01];
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PCAM, PRKMEANS, PRPLOTSOM, PREX_SOM

% Copyright: D. Tax, davidt@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

if nargin <5 | isempty(h)
	h = [0.6 0.2 0.01];
end
if nargin <4 | isempty(eta)
	eta = [0.5 0.3 0.1];
end
if nargin <3 | isempty(nrruns)
	nrruns = [20 40 40];
end
if nargin < 2| isempty(k),
	k = [5 5];
end
if nargin < 1 | isempty(a) 
    W = prmapping(mfilename,{k,nrruns,eta,h});
    W = setname(W,'Self-organising Map');
    return
end

if ~ismapping(k)           %training

	a = testdatasize(a);
	x = +a;     % use all the data!
	[nrx,dim] = size(x);
	% Set up the map size:
	if length(k)<2 % add a dummy dimension
		k = [k; 1];
	end

	% Some magic parameters:
	% The training consist of several runs, each differ in the width of
	% the neighborhood width, the learning rate and the number of SOM
	% updates:
	Dmax = distm(x,x); Dmax = max(Dmax(:));
	h = h*sqrt(Dmax);
	% (this should probably be something which can be set by the user)

	% Initialize the map
	w = 0.2*rand(k(1)*k(2),dim)+repmat(+mean(x),k(1)*k(2),1);

	% Run the different types of updates:
	for r1=1:length(nrruns)
		
		% Set up the neighbourhood:
		W = zeros(k(1),k(2));
		for i=1:k(1)
			for j=1:k(2)
				W(i,j) = exp(-((i-1)^2+(j-1)^2)/h(r1));
			end
		end

		% Do the real SOM training:
		for r2 = 1:nrruns(r1)

			%fprintf('%d/%d\t',r2,nrruns(r1));
			% Go randomly along the objects:
			I = randperm(nrx);
			for i=1:nrx
				% Find the winner:
				[mD,mI] = min(distm(w,x(I(i),:)),[],1);
				wx = rem(mI-1,k(1))+1;
				wy = floor((mI-1)/k(1))+1;
				% Update all weights:
				for j=1:k(1)*k(2)
					jx = rem(j-1,k(1))+1;
					jy = floor((j-1)/k(1))+1;
					w(j,:) = w(j,:) + eta(r1)*W(abs(jx-wx)+1,abs(jy-wy)+1)*...
						(x(I(i),:)-w(j,:));
				end
			end
		end  
	end

	clear W;
	% And save all useful data:
	V.k = k;  %(only used for plotting later...)
	V.neurons = w;
	W = prmapping(mfilename,'trained',V,(1:prod(k))',dim,prod(k));
	W = setname(W,'Self-organising Map');
else
    
	W = getdata(k); %unpack
	m = size(a,1); 
	% compute the distance to the nearest neuron in the map:
	%[mD,out] = min(distm(+a,W.w),[],2);
	D = distm(+a,W.neurons); % only compute the distances
	
	% Output a dataset with the negative distances:
	W = setdat(a,-D,k);
end











