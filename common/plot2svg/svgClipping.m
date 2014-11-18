function svgClipping(s, path)
% Adds a clipping path
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgClipping(s, path)
% Parameters:
%   s : Array of plot object handles
%   path : Clipping path nx3 or nx2.
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    userdata.svg.ClippingPath = path;
    set(s(i),'UserData', userdata);
end
