% PRIMPORT import the old-format prtools datasets
%
%   OUT = PRIMPORT(A)
%
% INPUT
%    A     The Structure to be converted.
%
% OUTPUT
%    OUT   The imported dataset
%
% DESCRIPTION 
% This routine converts old prtools datasets into the new prtools 4.x
% format. Structure A is tested for existence of all the fields forming
% particular dataset format. Prtools 3.x and 4.x formats are supported.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASET

% $Id: primport.m,v 1.4 2007/08/24 13:49:59 davidt Exp $

function out = primport(a)
	
	% we try to convert only structures
	if isstruct(a)

		% prtools 3.x format
		if isfield(struct(a),'d') & isfield(struct(a),'l') & ...
			    isfield(struct(a),'p') & isfield(struct(a),'f') & ...
			    isfield(struct(a),'ll') & isfield(struct(a),'c') & ...
			    isfield(struct(a), 's')
			
			prwarning(2,'prtools 3.x dataset converted to 4.x');
			
			% lab list could be a cell-array
			if iscell(a.ll)
				lablist=a.ll{1};
				nlab=a.l;
			else
				[nlab,lablist] = renumlab(a.l);
			end
			
			out = prdataset(a.d,nlab,'featlab',a.f,'prior',a.p,'lablist',lablist);
			
			if a.c<0
				% objects are pixels
				xx=-a.c;
				yy=size(a.d,1)/xx;
				out.objsize=[xx yy];
			end
			if a.c>0
				% objects are pictures
				xx=a.c;
				yy=size(a.d,2)/xx;
				out.featsize=[xx yy];
			end
			
		% prtools 4 format
		elseif isfield(struct(a),'data') & isfield(struct(a),'lablist') & ...
			    isfield(struct(a),'nlab') & isfield(struct(a),'labtype') & ...
			    isfield(struct(a),'targets') & isfield(struct(a), 'featlab') & ...
			    isfield(struct(a),'prior') & ...
			    isfield(struct(a),'objsize') & isfield(struct(a),'featsize') & ...
			    isfield(struct(a),'ident') & isfield(struct(a),'version') & ...
			    isfield(struct(a),'name') & isfield(struct(a),'user')
		
			%DXD: The previous line was commented out
			%     Uncommenting made it work:
			% out = prdataset(a.data,a.nlab,'lablist',a.lablist);
			out = prdataset(a.data);
			if ~isempty(a.nlab), out.nlab=a.nlab; end	
			if ~isempty(a.lablist), out.lablist=a.lablist; end	
			if ~isempty(a.labtype), out.labtype=a.labtype; end
			if ~isempty(a.targets), out.targets=a.targets; end
			if ~isempty(a.featlab), out.featlab=a.featlab; end
			if ~isempty(a.prior), out.prior=a.prior; end
			if ~isempty(a.objsize), out.objsize=a.objsize; end
			if ~isempty(a.featsize), out.featsize=a.featsize; end
			if ~isempty(a.ident), out.ident=a.ident; end
			if ~isempty(a.version), out.version=a.version; end
			if ~isempty(a.name), out.name=a.name; end
			if ~isempty(a.user), out.user=a.user; end

			% later, the featdom field was included
			if isfield(struct(a),'featdom') & ~isempty(a.featdom)
				out.featdom = a.featdom;
			end
			
			% later, the cost field was included
			if isfield(struct(a),'cost') & ~isempty(a.cost)
				out.cost = a.cost;
			end

			prwarning(2,'prtools 4.x converted from structure');
		
		else
			out = [];
			prwarning(2,'unknown format! The structure cannot be converted into a dataset');
		end
	
	else % anything else then structure: just copy to the output
		out=a;
	end
