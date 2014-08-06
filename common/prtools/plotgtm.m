%PLOTGTM  Plot a trained GTM mapping in 1D, 2D or 3D
%
%   H = PLOTGTM (W)
%
% INPUT
%   W   Trained GTM mapping
%
% OUTPUT
%   H   Graphics handles
%
% DESCRIPTION
% Creates a plot of the GTM manifold in the original data space, but at
% most in 3D.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% GTM, SOM, PRPLOTSOM

% (c) Dick de Ridder, 2003
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function ret = plotgtm (w)

        
    global GRIDSIZE;
    gs = [ 100 10 5 ]+1;
    
    if (~ismapping(w)) | (~strcmp(get(w,'name'),'GTM'))
        error ('can only plot a GTM mapping');
    else
        data = getdata(w); K = data{1}; M = data{2};
        W = data{3}; sigma = data{4}; mapto = data{5};
    end;
    
	[d,D] = size(w); KK = prod(K); MM = prod(M);
    
	phi_mu    = makegrid(M,D);				% Basis function centers.
	phi_sigma = 2/(mean(M)-1);				% Basis function widths.

    % For plotting.

    if (D == 1), 
		Kplot = gs(1); sz = 1/(gs(1)-1);
		xplot = 0:sz:1; 
	elseif (D == 2)
		Kplot = gs(2)^2; sz = 1/(gs(2)-1);
		[xx,yy] = meshgrid(0:sz:1,0:sz:1);
		xplot = [reshape(xx,1,Kplot); reshape(yy,1,Kplot)];
	elseif (D >= 3)
		Kplot = gs(3)^3; sz = 1/(gs(3)-1);
		[xx,yy,zz] = meshgrid(0:sz:1,0:sz:1,0:sz:1);
		xplot = [reshape(xx,1,Kplot); reshape(yy,1,Kplot); reshape(zz,1,Kplot)];
        D = 3;
	end;
 
	% Pre-calculate Phiplot.
    for j = 1:Kplot
    	for i = 1:MM
    		Phiplot(i,j) = exp(-(xplot(1:D,j)-phi_mu(1:D,i))'*(xplot(1:D,j)-phi_mu(1:D,i))/(phi_sigma^2));
    	end;
    end;

    % Plot grid lines.
    
    yplot = W(1:d,:)*Phiplot;

    hold on; h = [];
    if (D == 1)
	    h = [h plot(yplot(1,:),yplot(2,:),'k-')];
	elseif (D == 2)
    	for k = 1:gs(2)
    	    h = [h plot(yplot(1,(k-1)*gs(2)+1:k*gs(2)),yplot(2,(k-1)*gs(2)+1:k*gs(2)),'k-')];
    	    h = [h plot(yplot(1,k:gs(2):end),yplot(2,k:gs(2):end),'k-')];
    	end;
    elseif (D >= 3)
    	for k = 1:gs(3)
    	    for l = 1:gs(3)
        	    h = [h plot3(yplot(1,(k-1)*gs(3)^2+(l-1)*gs(3)+1:(k-1)*gs(3)^2+l*gs(3)),...
        	                 yplot(2,(k-1)*gs(3)^2+(l-1)*gs(3)+1:(k-1)*gs(3)^2+l*gs(3)),...
        	                 yplot(3,(k-1)*gs(3)^2+(l-1)*gs(3)+1:(k-1)*gs(3)^2+l*gs(3)),'k-')];
        	    h = [h plot3(yplot(1,(k-1)*gs(3)^2+l:gs(3):k*gs(3)^2),...
        	                 yplot(2,(k-1)*gs(3)^2+l:gs(3):k*gs(3)^2),...
        	                 yplot(3,(k-1)*gs(3)^2+l:gs(3):k*gs(3)^2),'k-')];
        	    h = [h plot3(yplot(1,(k-1)*gs(3)+l:gs(3)^2:end),...
        	                 yplot(2,(k-1)*gs(3)+l:gs(3)^2:end),...
        	                 yplot(3,(k-1)*gs(3)+l:gs(3)^2:end),'k-')];
            end;
        end;
 	    view(3);
	end;
    set (h,'LineWidth',2);
	hold off;
    
    if (nargout > 0), ret = h; end;
    
return
	
% GRID = MAKEGRID (K,D)
%
% Create a KK = prod(K)-dimensional grid with dimensions K(1), K(2), ... of 
% D-dimensional uniformly spaced grid points on [0,1]^prod(K), and store it 
% as a D x KK matrix X.

function grid = makegrid (K,D)

	KK = prod(K);

	for h = 1:D
		xx{h} = 0:(1/(K(h)-1)):1;										% Support point means
	end;

	% Do that voodoo / that you do / so well...

  if (D==1)
  	xm = xx;
  else
  	cmd = '[';
    for h = 1:D-1, cmd = sprintf ('%sxm{%d}, ', cmd, h); end; 
    cmd = sprintf ('%sxm{%d}] = ndgrid(', cmd, D);
    for h = 1:D-1, cmd = sprintf ('%sxx{%d}, ', cmd, h); end; 
    cmd = sprintf ('%sxx{%d});', cmd, D); eval(cmd);
  end;
    
	cmd = 'mm = zeros(D, ';
	for h = 1:D-1, cmd = sprintf ('%s%d, ', cmd, K(h)); end; 
	cmd = sprintf ('%s%d);', cmd, K(D)); eval (cmd);

	for h = 1:D
		cmd = sprintf ('mm(%d,', h);
		for g = 1:D-1, cmd = sprintf ('%s:,', cmd); end; 
    cmd = sprintf ('%s:) = xm{%d};', cmd, h); eval (cmd);
	end;

	grid = reshape(mm,D,KK);

return



