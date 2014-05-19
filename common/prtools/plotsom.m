function h = plotsom(W)
%PLOTSOM Plot the Self-Organizing Map in 2D
%
%    PLOTSOM(W)
%
% Plot the Self-Organizing Map W, trained by som.m. This is only
% possible if the map is 2D.
%
% SEE ALSO
% SOM

% Copyright: D.M.J. Tax, davidt@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% Maybe I should introduce the possibility to set the linewidth and
% markersize... 
% 
% Changes:
% DR1 - Dick de Ridder, 02-03-2006:
%       If the mapped is trained on 3D data, plot it like this.

if ~ismapping(W)
	error('I expect a SOM mapping!');
end
if ( ~strcmp(getmapping_file(W),'som') & ...
	 ~strcmp(getmapping_file(W),'som_dd') )
    error('I expect a SOM mapping!');
end
if size(W.data.neurons,2)~=2
	error('The SOM can only be plotted in 2D');
end

% Get the data:
W = +W;  w=W.neurons; k=W.k;
% Plot the bloody thing:
hold on;

% DR: Handle maps trained on 3D data.

if (size(w,2) == 3)
	% The 'horizontal' lines:
	for i=0:k(2)-1
			h=plot3(w(i*k(1)+(1:k(1)),1),w(i*k(1)+(1:k(1)),2),w(i*k(1)+(1:k(1)),3),'o-');
			set(h,'linewidth',2,'markersize',8);
	end
	I = reshape(1:k(1)*k(2),k(1),k(2))';
	% The 'vertical' lines:
	for i=0:k(1)-1
			h=plot3(w(I(i*k(2)+(1:k(2))),1),w(I(i*k(2)+(1:k(2))),2),w(I(i*k(2)+(1:k(2))),3),'o-');
			set(h,'linewidth',2,'markersize',8);
	end
	view(3);
else
	% The 'horizontal' lines:
	for i=0:k(2)-1
			h=plot(w(i*k(1)+(1:k(1)),1),w(i*k(1)+(1:k(1)),2),'o-');
			set(h,'linewidth',2,'markersize',8);
	end
	I = reshape(1:k(1)*k(2),k(1),k(2))';
	% The 'vertical' lines:
	for i=0:k(1)-1
			h=plot(w(I(i*k(2)+(1:k(2))),1),w(I(i*k(2)+(1:k(2))),2),'o-');
			set(h,'linewidth',2,'markersize',8);
	end
end;

% Return the handle only when it is required
if (nargout==0)
	clear h;
end

return
