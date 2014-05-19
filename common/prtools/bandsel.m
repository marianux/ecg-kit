%BANDSEL Selection of bands from object images
%
%   B = BANDSEL(A,J)
%   W = BANDSEL([],J)
%   B = A*BANDSEL([],J)
%
% INPUT
%   A    Dataset or datafile with multi-band object images
%   J    Indices of bands to be selected
%
% OUTPUT
%   W    Mapping performing the band selection
%   B    Dataset with selected bands (ordered according to J)
%
% DESCRIPTION
% If the objects in a dataset or datafile A are multi-band images, e.g. RGB
% images, or the result of IM_PATCH, then the featsize of A is [M,N,L] for
% for L bands of an M x N images. This routine makes a selection J out of
% L. The routine BAND2OBJ may be used to organize the bands vertically
% as separate objects. However, BANDSEL nor BAND2OBJ can be applied to 
% datafiles for which already a bandselection has been defined by BANDSEL.
%
%
% SEE ALSO
% DATASETS, DATAFILES, IM2OBJ, IM_PATCH, BAND2OBJ

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = bandsel(a,J)

mapname = 'Band Selection';
if nargin < 2, J = 1; end

if nargin < 1 | isempty(a)

	w = prmapping(mfilename,'fixed',{J});
	w = setname(w,mapname);

elseif isdatafile(a) % just store administration

	if size(J,1) ~= 1
		if size(J,1) ~= size(a,1)
			error('Matrix with band indices does not match number of objects')
		end
	else
		J = repmat(J,size(a,1),1);
	end
	
	%determine number of bands to be selected
	if iscell(J)
		n = size(J{1},2);
	else
		n = size(J,2);
	end
	
	%store bandselection as a mapping to be executed 
	%during dataset conversion
	v = prmapping(mfilename,'fixed',{[]},[],0,0);

	% make the mappings identical, i.e. independent of v.data
	% and store the data in the ident field under 'bandsel'.
	% This will enable vertical concatenation of datafiles

	J0 = getident(a,'bandsel');
	J1 = zeros(size(J));
	if ~isempty(J0)
		% we already have a bandselection set; change it
		if size(J,1) == 1
			J1 = J0(:,J);
		else
			for j=1:size(J0,1)
				J1(j,:) = J0(j,J(j,:));
			end
		end
	else
		J1 = J;
	end
	
	%correct bandnames in ident
	a = setident(a,J1,'bandsel');
	bandnames = getident(a,'bandnames');
	
	if ~isempty(bandnames)
		for j=1:size(a,1)
			bandnames{j} = bandnames{j}(J(j,:),:);
		end
		a = setident(a,bandnames,'bandnames');
	end
	v = setdata(v,[]);
	w = addpostproc(a,v);

elseif isdataset(a) % execute

	m = size(a,1);
	if isempty(J)  % J is stored in a.ident.bandsel
		J = getident(a,'bandsel');
		% we assume that new bandnames are already set,
		% simulataneously with bandsel
		if isempty(J) % no bandselection defined
			w = a;
			return
		end
	else
		if size(J,1) == 1
			J = repmat(J(:)',m,1);
		else
			if size(J,1) ~= m
				error('Wrong size of band selection array')
			end
		end
		bandnames = getident(a,'bandnames');
		if ~isempty(bandnames)
			for j=1:m
				bandnames{j} = bandnames{j}(J,:);
			end
			a = setident(a,bandnames,'bandnames');
		end
	end
	isobjim(a);
	fsize = getfeatsize(a);
	if length(fsize) < 3
		error('No image bands defined for dataset')
	end
	if any(J(:) > fsize(3))
		id = getident(a);
		id = id(find(J>fsize(3)));
		error(['Wrong index for image bands. Object ident: ' int2str(id(1))] )
	end
	k = prod(fsize(1:2));                           % size of a band
	L = repmat((J(:,1)-1)*k,1,k)+repmat([1:k],m,1); % indices of band_1 per object
	if size(J,2) > 1                                % concatenate in case of multiple band selection
		for j=2:size(J,2)
			LL = repmat((J(:,j)-1)*k,1,k)+repmat([1:k],m,1); % indices of band_j per object
			L = [L LL];                                % indices for all bands to be selected
		end
	end
	adata = getdata(a);                            % Let us do the selection on the data
	bdata = zeros(m,k*size(J,2));
	for i=1:m                                      % can this be done faster?
		bdata(i,:) = adata(i,L(i,:));
	end
	b = setdat(a,bdata);                           % store the data in a dataset with all information of a
	b = setident(b,[],'bandsel');                  % bandselection done, avoid second in case of stacked selections
	w = setfeatsize(b,[fsize(1:2) size(J,2)]); 

else

	error('Illegal command')

end
