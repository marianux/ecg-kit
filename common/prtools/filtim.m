%FILTIM Mapping to filter multiband image objects in datasets and datafiles
%
%    B = FILTIM(A,FILTER_COMMAND,{PAR1,PAR2,....},SIZE)
%    B = A*FILTIM([],FILTER_COMMAND,{PAR1,PAR2,....},SIZE)
%
% INPUT
%    A               Dataset or datafile with multi-band image objects
%    FILTER_COMMAND  String with function name
%    {PAR1, ...  }   Cell array with optional parameters to FILTER_COMMAND
%    SIZE            Output size of the mapping (default: input size)
%
% OUTPUT
%    B               Dataset containing multi-band images processed by 
%                    FILTER_COMMAND, band by band. 
%
% DESCRIPTION
% For each band of each object stored in A a filter operation is performed as
%
%    OBJECT_OUT = FILTER_COMMAND(OBJECT_IN,PAR1,PAR2,....)
%
% The results are collected and stored in B. In case A (and thereby B) is
% a datafile, execution is postponed until conversion into a dataset, or a
% call to SAVEDATAFILE.
%
% The difference between FILTIM and the similar command FILTM is that
% FILTIM is aware of the band structure of the objects. As FILTIM treats
% the bands separately it cannot be used for commands that change the number
% of bands (like RGB2GRAY) or need to access them all.
%
% EXAMPLE
% a = delft_images; b = a(120 121 131 230)*col2gray
% e = b*filtim([],'fft2')*filtim([],'abs')*filtim([],'fftshift');
% figure; show(e); figure; show((1+e)*filtim([],'log')); 
%
% SEE ALSO
% DATASETS, DATAFILES, IM2OBJ, DATA2IM, IM2FEAT, DATGAUSS, DATFILT, FILTM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = filtim(a,command,pars,outsize)

		if nargin < 4, outsize = []; end
	if nargin < 3, pars = {}; end
	if nargin < 2
		error('No command given')
	end
	if ~iscell(pars), pars = {pars}; end
	
	mapname = 'Dataset/file band filtering';
	
	if isempty(a)                    % no data, so just mapping definition
		b = prmapping(mfilename,'fixed',{command,pars});
		if ~isempty(outsize)
			b = setsize_out(b,outsize);
		end
		b = setname(b,mapname);

	elseif isdatafile(a)                 
		
		% for datafiles filters will be stored
		isobjim(a);      
    if isempty(getpostproc(a))     % as preprocessing (if no postproc defined)
      b = addpreproc(a,mfilename,{command pars},outsize);
    else                           % or as mapping as postprocessing
		  v = prmapping(mfilename,'fixed',{command,pars});
		  if ~isempty(outsize)
			  v = setsize_out(v,outsize);
		  end
		  v = setname(v,mapname);
		  b = addpostproc(a,v);
		end
    
	elseif isdataset(a)
   
		% convert to image and process
		isobjim(a);
    m = size(a,1);
    d = data2im(a);
		out = feval(mfilename,d,command,pars); % effectively jumps to double (below)
    out = double(out);
		fsize = size(out);		
		if m > 1
			fsize = fsize(1:3);
			if fsize(3) == 1
				fsize = fsize(1:2);
			end
		end
    % store processed images in dataset
		b = setdata(a,im2obj(out,fsize));
		b = setfeatsize(b,fsize);
		
	else % double
		
		% make imsize 4D: horz*vert*bands*objects
    imsize = size(a);
		if length(imsize) == 2
			imsize = [imsize 1 1];
		elseif length(imsize) == 3
			imsize = [imsize 1];
		end
    
		% find size of first image
		
		
		first = execute(command,getsubim(a,1,1),pars);		
	  %first = feval(command,getsubim(a,1,1),pars{:});
		if all(imsize(3:4) == 1) % preserve DipLib format for the time being
			b = first;
			return
		end
    first = double(first); % for Dip_Image users
		outsize = size(first);
		if length(outsize) == 3 % single subimage generates multiple bands !!
			b = repmat(first(:,:,1),[1 1 imsize(3)*outsize(3) imsize(4)]);
		else
			b = repmat(first,[1 1 imsize(3:4)]);
		end
    % process all other images
		[nn,s,count] = prwaitbarinit('Filtering %i images',imsize(4));
		for j = 1:imsize(4)
			for i = 1:imsize(3)
				ima = double(execute(command,getsubim(a,i,j),pars));
				%ima = double(feval(command,getsubim(a,i,j),pars{:}));
				if (any(outsize ~= size(ima)))  % check size
					error('All image sizes should be the same')
		  	end
				if length(outsize) == 2 % simple case: image_in --> image_out
					b(:,:,i,j) = ima;
				else                    % advanced: image_in --> bands out
					b(:,:,i:imsize(3):end,j) = ima;
				end
			end
			count = prwaitbarnext(imsize(4),'Filtering %i images: ',count);
		end
  end
  
return

function asub = getsubim(a,i,j)

% needed as Dip_Image cannot handle 4d subsript if data is 2d or 3d

n = ndims(a);

if n == 2
	if i ~= 1 | j ~= 1
		error('Illegal image requested')
	end
	asub = a;
elseif n == 3
	if j ~= 1
		error('Illegal image requested')
	end
	asub = a(:,:,i);
elseif n == 4
	asub = a(:,:,i,j);
else
	error('Incorrect image supplied')
end

function out = execute(command,a,pars)

	exist_command = exist(command);
	if isstr(command) & any([2 3 5 6] == exist_command)
		if isempty(pars) | numel(pars) == 1 & isempty(pars{:})
			out = feval(command,a);
		else
			out = feval(command,a,pars{:});
		end
	elseif ismapping(command)
		out = prmap(a,command);
	else
		error([command ': Filter command not found'])
	end
 
return
		

 
