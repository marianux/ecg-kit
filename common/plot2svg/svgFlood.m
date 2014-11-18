function svgFlood(s, color, opacity, result)
% Adds a feFlood SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgFlood(s, color, opacity, result)
% Parameters:
%   s : Array of plot object handles
%   color  : A 1x3 color vector [0...1[
%   opacity : Opacity value [0...1]
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feFlood';
    userdata.svg.Filter(next).Subfilter.Color = color;
    userdata.svg.Filter(next).Subfilter.Opacity = opacity;
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end
